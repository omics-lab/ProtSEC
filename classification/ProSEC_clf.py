import numpy as np
from embedder import ProteinEmbedder


import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
from sklearn.metrics import accuracy_score, classification_report
import matplotlib.pyplot as plt


# Check for device availability
if torch.backends.mps.is_available():
    device = torch.device("mps")
    print("Using MPS (Apple Silicon)")
elif torch.cuda.is_available():
    device = torch.device("cuda")
    print("Using CUDA")
else:
    device = torch.device("cpu")
    print("Using CPU")


class AttentionModule(nn.Module):
    def __init__(self, hidden_dim):
        super(AttentionModule, self).__init__()
        self.attention = nn.Linear(hidden_dim, 1)
        
    def forward(self, lstm_output):
        # lstm_output shape: (batch_size, seq_len, hidden_dim)
        attention_weights = torch.softmax(self.attention(lstm_output), dim=1)
        # attention_weights shape: (batch_size, seq_len, 1)
        
        # Apply attention weights
        attended_output = torch.sum(attention_weights * lstm_output, dim=1)
        # attended_output shape: (batch_size, hidden_dim)
        
        return attended_output, attention_weights


class CNNLSTMClassifier(nn.Module):
    def __init__(self, input_size=2000, embedding_dim=20, num_classes=10):
        super(CNNLSTMClassifier, self).__init__()
        
        self.input_size = input_size
        self.embedding_dim = embedding_dim
        
        # Input embedding layer to transform 2000 features to (seq_len, embedding_dim)
        # We'll treat the 2000 features as a sequence of length 2000 with 1 feature each
        # Then embed to embedding_dim
        self.input_embedding = nn.Linear(1, embedding_dim)
        
        # Convolutional layers
        self.conv1 = nn.Conv1d(embedding_dim, 10, kernel_size=3, padding=1)
        self.conv2 = nn.Conv1d(embedding_dim, 10, kernel_size=5, padding=2)
        self.final_conv = nn.Conv1d(20, 20, kernel_size=3, padding=1)
        
        # LSTM layers
        self.forward_lstm = nn.LSTM(20, 15, batch_first=True)
        self.backward_lstm = nn.LSTM(20, 15, batch_first=True)
        
        # Attention mechanism
        self.attention = AttentionModule(30)  # 15 + 15 from bidirectional LSTM
        
        # Dense layers
        self.dense = nn.Linear(30, 30)
        self.dropout = nn.Dropout(0.3)
        
        # Output layer
        self.output = nn.Linear(30, num_classes)
        
    def forward(self, x):
        batch_size = x.size(0)
        
        # Reshape input: (batch_size, 2000) -> (batch_size, 2000, 1)
        x = x.unsqueeze(-1)
        
        # Apply input embedding: (batch_size, 2000, 1) -> (batch_size, 2000, 20)
        x = self.input_embedding(x)
        
        # Transpose for conv1d: (batch_size, 2000, 20) -> (batch_size, 20, 2000)
        x = x.transpose(1, 2)
        
        # Apply convolutional layers
        conv1_out = F.relu(self.conv1(x))  # (batch_size, 10, 2000)
        conv2_out = F.relu(self.conv2(x))  # (batch_size, 10, 2000)
        
        # Concatenate conv outputs
        conv_concat = torch.cat([conv1_out, conv2_out], dim=1)  # (batch_size, 20, 2000)
        
        # Final conv layer
        conv_final = F.relu(self.final_conv(conv_concat))  # (batch_size, 20, 2000)
        
        # Transpose back for LSTM: (batch_size, 20, 2000) -> (batch_size, 2000, 20)
        x = conv_final.transpose(1, 2)
        
        # Forward LSTM
        forward_out, _ = self.forward_lstm(x)  # (batch_size, 2000, 15)
        
        # Backward LSTM (reverse the sequence)
        x_reversed = torch.flip(x, dims=[1])
        backward_out, _ = self.backward_lstm(x_reversed)  # (batch_size, 2000, 15)
        backward_out = torch.flip(backward_out, dims=[1])  # Flip back
        
        # Concatenate LSTM outputs
        lstm_concat = torch.cat([forward_out, backward_out], dim=2)  # (batch_size, 2000, 30)
        
        # Apply attention
        attended_output, attention_weights = self.attention(lstm_concat)  # (batch_size, 30)
        
        # Dense layer
        dense_out = F.relu(self.dense(attended_output))
        dense_out = self.dropout(dense_out)
        
        # Output layer
        output = self.output(dense_out)
        
        return output

def train_model(model, train_loader, val_loader, num_epochs=50, learning_rate=0.001):
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5)
    
    train_losses = []
    val_losses = []
    train_accuracies = []
    val_accuracies = []
    
    best_val_acc = 0.0
    best_model_state = None
    
    for epoch in range(num_epochs):
        # Training phase
        model.train()
        train_loss = 0.0
        train_correct = 0
        train_total = 0
        
        for batch_x, batch_y in train_loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            
            optimizer.zero_grad()
            outputs = model(batch_x)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
            _, predicted = torch.max(outputs.data, 1)
            train_total += batch_y.size(0)
            train_correct += (predicted == batch_y).sum().item()
        
        # Validation phase
        model.eval()
        val_loss = 0.0
        val_correct = 0
        val_total = 0
        
        with torch.no_grad():
            for batch_x, batch_y in val_loader:
                batch_x, batch_y = batch_x.to(device), batch_y.to(device)
                outputs = model(batch_x)
                loss = criterion(outputs, batch_y)
                
                val_loss += loss.item()
                _, predicted = torch.max(outputs.data, 1)
                val_total += batch_y.size(0)
                val_correct += (predicted == batch_y).sum().item()
        
        # Calculate metrics
        train_acc = 100 * train_correct / train_total
        val_acc = 100 * val_correct / val_total
        avg_train_loss = train_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)
        
        train_losses.append(avg_train_loss)
        val_losses.append(avg_val_loss)
        train_accuracies.append(train_acc)
        val_accuracies.append(val_acc)
        
        # Save best model
        if val_acc > best_val_acc:
            best_val_acc = val_acc
            best_model_state = model.state_dict().copy()
        
        scheduler.step(avg_val_loss)
        
        if epoch % 5 == 0:
            print(f'Epoch [{epoch+1}/{num_epochs}]')
            print(f'Train Loss: {avg_train_loss:.4f}, Train Acc: {train_acc:.2f}%')
            print(f'Val Loss: {avg_val_loss:.4f}, Val Acc: {val_acc:.2f}%')
            print('-' * 50)
    
    # Load best model
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
    
    return model, train_losses, val_losses, train_accuracies, val_accuracies

def evaluate_model(model, test_loader):
    model.eval()
    all_predictions = []
    all_labels = []
    
    with torch.no_grad():
        for batch_x, batch_y in test_loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            outputs = model(batch_x)
            _, predicted = torch.max(outputs.data, 1)
            
            all_predictions.extend(predicted.cpu().numpy())
            all_labels.extend(batch_y.cpu().numpy())
    
    accuracy = accuracy_score(all_labels, all_predictions)
    report = classification_report(all_labels, all_predictions)
    
    return accuracy, report, all_predictions, all_labels

def plot_training_history(train_losses, val_losses, train_accuracies, val_accuracies):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # Plot losses
    ax1.plot(train_losses, label='Training Loss')
    ax1.plot(val_losses, label='Validation Loss')
    ax1.set_title('Model Loss')
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('Loss')
    ax1.legend()
    
    # Plot accuracies
    ax2.plot(train_accuracies, label='Training Accuracy')
    ax2.plot(val_accuracies, label='Validation Accuracy')
    ax2.set_title('Model Accuracy')
    ax2.set_xlabel('Epoch')
    ax2.set_ylabel('Accuracy (%)')
    ax2.legend()
    
    plt.tight_layout()
    plt.show()




def cut_protein_sequence_center(protein_sequence, n):
    """
    Cut protein sequence from the center to achieve desired length n.
    
    If the sequence length is > n, it will remove the central portion
    and keep the N-terminal and C-terminal ends.
    
    Args:
        protein_sequence (str): Input protein sequence
        n (int): Desired output length
    
    Returns:
        str: Protein sequence of length n (if original length > n)
              or original sequence (if original length <= n)
    
    Examples:
        >>> cut_protein_sequence_center("ABCDEFGHIJ", 6)
        'ABCHIJ'
        >>> cut_protein_sequence_center("ABCDEFGHIJ", 4)
        'ABHIJ'
        >>> cut_protein_sequence_center("ABC", 6)
        'ABC'
    """
    # If sequence is already shorter than or equal to n, return as is
    if len(protein_sequence) <= n:
        return protein_sequence
    
    # Calculate how many residues to keep from each end
    total_to_keep = n
    
    # Split the kept residues equally between N-terminal and C-terminal
    # If odd number to keep, keep one more from the N-terminal
    n_terminal_keep = (total_to_keep) // 2
    c_terminal_keep = total_to_keep - n_terminal_keep
    
    # Extract N-terminal and C-terminal portions
    n_terminal = protein_sequence[:n_terminal_keep]
    c_terminal = protein_sequence[-c_terminal_keep:] if c_terminal_keep > 0 else ""
    
    # Combine N-terminal and C-terminal portions
    return n_terminal + c_terminal

def load_data_from_npz(file_path):
    """
    Load data from an NPZ file.
    
    Args:
        file_path (str): NPZ file to load
        
    Returns:
        tuple: (X_train, Y_train, X_test, Y_test)
    """
    try:
        data = np.load(file_path, allow_pickle=True)
        X_train = data['X_train']
        Y_train = data['Y_train']
        X_test = data['X_test']
        Y_test = data['Y_test']
        
        print(f"Data successfully loaded from '{file_path}'")
        return X_train, Y_train, X_test, Y_test
        
    except Exception as e:
        print(f"Error loading data from NPZ file: {e}")
        return None, None, None, None




if __name__ == "__main__":
    clf_data_path = 'classification/Deeploc/clf_data.npz'
    X_train, Y_train, X_test, Y_test = load_data_from_npz(clf_data_path)
    if X_train is not None and Y_train is not None:
        print(f"Training data shape: {X_train.shape}, Labels shape: {Y_train.shape}")
        print(f"Test data shape: {X_test.shape}, Labels shape: {Y_test.shape}")
    else:
        print("Failed to load data from NPZ file.")
    
    # Initialize the protein embedder
    embedder = ProteinEmbedder(dim_reduct='MDS', dist_func='SMS', n=1000)

    X_train_embedded = np.array([embedder.encode(cut_protein_sequence_center(seq, 1000)) for seq in X_train])
    X_test_embedded = np.array([embedder.encode(cut_protein_sequence_center(seq, 1000)) for seq in X_test])

    # Convert to numpy arrays
    X_train_embedded = np.concatenate([X_train_embedded.real, X_train_embedded.imag], axis=1)
    X_test_embedded = np.concatenate([X_test_embedded.real, X_test_embedded.imag], axis=1)

    print(f"Embedded training data shape: {X_train_embedded.shape}")
    print(f"Embedded test data shape: {X_test_embedded.shape}")

    unique_labels = list(set(Y_train))

    Y_train_encoded = np.array([unique_labels.index(label) for label in Y_train])
    Y_test_encoded = np.array([unique_labels.index(label) for label in Y_test])

    

    X_train_tensor = torch.FloatTensor(X_train_embedded)
    X_test_tensor = torch.FloatTensor(X_test_embedded)
    Y_train_tensor = torch.LongTensor(Y_train_encoded)
    Y_test_tensor = torch.LongTensor(Y_test_encoded)

    # Create datasets and dataloaders
    train_dataset = TensorDataset(X_train_tensor, Y_train_tensor)
    test_dataset = TensorDataset(X_test_tensor, Y_test_tensor)

    # Split training data for validation
    val_size = int(0.1 * len(train_dataset))
    train_size = len(train_dataset) - val_size
    train_subset, val_subset = torch.utils.data.random_split(train_dataset, [train_size, val_size])

    train_loader = DataLoader(train_subset, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_subset, batch_size=64, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=64, shuffle=False)

    # Initialize model
    model = CNNLSTMClassifier(input_size=2000, embedding_dim=16, num_classes=10).to(device)

    print(f"Model parameters: {sum(p.numel() for p in model.parameters() if p.requires_grad)}")

    # Train the model
    model, train_losses, val_losses, train_accs, val_accs = train_model(
        model, train_loader, val_loader, num_epochs=150, learning_rate=0.01
    )

    # Evaluate on test set
    test_accuracy, test_report, predictions, true_labels = evaluate_model(model, test_loader)

    print(f"Test Accuracy: {test_accuracy:.4f}")
    print("\nClassification Report:")
    print(test_report)

