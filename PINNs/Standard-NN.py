from PIL import Image
import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt

# Helper function for saving GIFs
def save_gif_PIL(outfile, files, fps=5, loop=0):
    imgs = [Image.open(file) for file in files]
    imgs[0].save(fp=outfile, format='GIF', append_images=imgs[1:], save_all=True, duration=int(1000/fps), loop=loop)

# Analytical solution of the underdamped harmonic oscillator
def oscillator(d, w0, x):
    assert d < w0
    w = np.sqrt(w0**2 - d**2)
    phi = np.arctan(-d / w)
    A = 1 / (2 * np.cos(phi))
    cos = torch.cos(phi + w * x)
    exp = torch.exp(-d * x)
    y = exp * 2 * A * cos
    return y

# Define a fully connected neural network
class FCN(nn.Module):
    def __init__(self, N_INPUT, N_OUTPUT, N_HIDDEN, N_LAYERS):
        super().__init__()
        activation = nn.Tanh
        self.fcs = nn.Sequential(nn.Linear(N_INPUT, N_HIDDEN), activation())
        self.fch = nn.Sequential(*[nn.Sequential(nn.Linear(N_HIDDEN, N_HIDDEN), activation()) for _ in range(N_LAYERS - 1)])
        self.fce = nn.Linear(N_HIDDEN, N_OUTPUT)

    def forward(self, x):
        x = self.fcs(x)
        x = self.fch(x)
        x = self.fce(x)
        return x

# Generate training data
d, w0 = 2, 20
x = torch.linspace(0, 1, 500).view(-1, 1)  # Full domain
y = oscillator(d, w0, x).view(-1, 1)
x_data = x[0:200:20]  # Training points
y_data = y[0:200:20]

# Plot exact solution and training points
plt.figure()
plt.plot(x, y, label="Exact solution")
plt.scatter(x_data, y_data, color="tab:orange", label="Training data")
plt.legend()
plt.show()

# Helper function to plot training progress
def plot_result(x, y, x_data, y_data, yh, xp=None, i=0):
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, color="grey", linewidth=2, alpha=0.8, label="Exact solution")
    plt.plot(x, yh, color="tab:blue", linewidth=4, alpha=0.8, label="Neural network prediction")
    plt.scatter(x_data, y_data, s=60, color="tab:orange", alpha=0.4, label="Training data")
    if xp is not None:
        plt.scatter(xp, -0 * torch.ones_like(xp), s=60, color="tab:green", alpha=0.4, label="Physics loss training locations")
    l = plt.legend(loc=(1.01, 0.34), frameon=False, fontsize="large")
    plt.setp(l.get_texts(), color="k")
    plt.xlim(-0.05, 1.05)
    plt.ylim(-1.1, 1.1)
    plt.text(1.065, 0.7, "Training step: %i" % (i + 1), fontsize="xx-large", color="k")
    plt.axis("off")

# Train the neural network
torch.manual_seed(123)
model = FCN(1, 1, 32, 3)  # Define model
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)  # Define optimizer

files = []
for i in range(1000):
    optimizer.zero_grad()
    yh = model(x_data)  # Predict on training data
    loss = torch.mean((yh - y_data) ** 2)  # Mean squared error loss
    loss.backward()
    optimizer.step()

    # Plot and save training progress every 10 steps
    if (i + 1) % 10 == 0:
        yh = model(x).detach()
        plot_result(x, y, x_data, y_data, yh, i=i)
        file = f"plots/nn_{i + 1:08d}.png"
        plt.savefig(file, bbox_inches='tight', pad_inches=0.1, dpi=100, facecolor="white")
        files.append(file)
        if (i + 1) % 500 == 0:
            plt.show()
        else:
            plt.close("all")

# Save training progress as a GIF
save_gif_PIL("nn.gif", files, fps=20, loop=0)

