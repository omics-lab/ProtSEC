## Recommendations for selecting the embedding dimension

The embedding vector length (d) determines the output length of the FFT. Appropriate selection of d is essential for balancing computational efficiency with information preservation. The relationship between d and input sequence length follows three cases:

- When d < input length, sequences are truncated from the end.
- When d > input length, sequences are zero-padded at the end.
- When d is unspecified, a default value of d = 1024 is applied.

### For datasets with substantial sequence length variation

For datasets exhibiting substantial variation in sequence lengths, we recommend setting d to the third quartile (Q3) of the length distribution. To optimize FFT computational performance, d can be rounded up to the nearest power of two (≥ Q3), as FFT algorithms achieve maximum efficiency at power-of-two lengths.

The implications of this choice include:

- Sequences exceeding d will be truncated, resulting in minor information loss.
- Sequences shorter than d will be zero-padded.

This strategy represents a practical trade-off between computational speed and the marginal information loss incurred through truncation.

### For datasets with uniformly short sequences

For datasets characterized by uniformly short sequences, d should be set to the maximum observed sequence length. This eliminates truncation entirely, ensuring complete information preservation. Rounding d up to the nearest power of two further enhances FFT efficiency.

The implications of this choice include:

- All sequences will be zero-padded to length d.
- FFT computation remains efficient while preserving all original sequence information.
