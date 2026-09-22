import numpy as np
import jax
import jax.numpy as jnp

from utility import aa_to_complex

from Bio import SeqIO


# JAX defaults to 32-bit. True matches the original complex128 output;
# set to False for complex64 (half the memory, noticeably faster, esp. on GPU).
# USE_X64 = True
# jax.config.update("jax_enable_x64", USE_X64)

PAD = 0  # byte code used for padding (NUL never appears in a sequence) -> maps to 0j


class ProteinEmbedder:

    def __init__(self, dim_reduct='MDS', dist_func='ASMP', n=256):
        """
        Args:
            dim_reduct: Dimensionality-reduction method used for the mapping ('t-SNE', 'MDS', 'UMAP')
            dist_func:  Distance function used for the mapping
            n:          Fixed vector length; shorter sequences are zero-padded, longer ones truncated
        """
        # Copy so the shared dict in `utility` is not modified
        mapping = dict(aa_to_complex[dim_reduct][dist_func])
        mapping['X'] = sum(mapping.values()) / len(mapping)
        self.aa_to_complex = mapping
        self.n = n

        # 256-entry lookup table indexed by byte value.
        # Every unknown character defaults to 'X'; upper and lower case both handled.
        lut = np.full(256, mapping['X'], dtype=np.complex64)
        for aa, z in mapping.items():
            lut[ord(aa.upper())] = z
            lut[ord(aa.lower())] = z
        lut[PAD] = 0j
        self.lut = jnp.asarray(lut)

        # note: self.lut is captured as a closed-over constant inside embed_one,
        # and _jitted_embed is compiled against that constant right here.
        # Treat self.lut as immutable after this point — if it were ever
        # reassigned post-construction, this compiled function would keep
        # using the stale values. Use set_lut() below if you need to change it.
        self._build_jitted_embed()

    def _build_jitted_embed(self):
        """(Re)compile the vmapped+jitted embedding pipeline against the
        current self.lut. Called from __init__, and again from set_lut()
        if the lookup table is ever replaced.
        """
        self._jitted_embed = jax.jit(jax.vmap(self.embed_one))

    def set_lut(self, new_lut):
        """Safely replace the lookup table and recompile the jitted pipeline,
        so embed_many() never silently uses a stale closed-over lut.
        """
        self.lut = jnp.asarray(new_lut, dtype=jnp.complex64)
        self._build_jitted_embed()

    def to_codes(self, sequence):
        """Part 1 (plain Python/NumPy): string -> fixed-length uint8 codes.
        Must run OUTSIDE jax.vmap/jit, since it does string slicing/encoding
        that JAX tracers can't handle.
        """
        codes = np.frombuffer(sequence[:self.n].encode('ascii', errors='replace'), dtype=np.uint8)
        return np.pad(codes, (0, self.n - len(codes)), constant_values=0)

    def embed_one(self, codes):
        return jnp.fft.fft(self.lut[codes])

    def embed_many(self, sequences):
        """Part 2 (JAX): the per-batch function, numeric only.

        `sequences` is a list/iterable of plain Python strings.

        note on recompilation: codes always has shape (batch, self.n), since
        each sequence is padded/truncated to a fixed length in to_codes().
        JAX/XLA will trigger a fresh trace+compile of _jitted_embed the first
        time a given `batch` size is seen, then reuse that compiled version
        for subsequent calls with the same batch size. If you call this
        repeatedly with varying numbers of sequences per call, you'll pay a
        recompilation cost each time a new batch size shows up. If that
        matters for your workload, call in fixed-size chunks (pad the last
        partial chunk) so only a handful of shapes are ever compiled.
        """
        # String -> code conversion happens in plain Python/NumPy, not under vmap.
        codes = np.stack([self.to_codes(s) for s in sequences])
        spectra = self._jitted_embed(codes)  # embed_one already does lut-lookup + fft
        return spectra


if __name__ == '__main__':
    import time
    embedder = ProteinEmbedder()
    # curl -O https://services.healthtech.dtu.dk/services/DeepLoc-1.0/deeploc_data.fasta
    
    fasta_path = './deeploc_data.fasta'

    sequences = [str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")] * 72
    print(f'Total sequences: {len(sequences)}')


    # for record in SeqIO.parse(fasta_path, "fasta"):
    #     print(record.id)

    t0 = time.perf_counter()

    spectra = embedder.embed_many(sequences)
    t1 = time.perf_counter()
    print(f"Time taken: {t1 - t0:.3f} seconds")
