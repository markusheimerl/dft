import torch
import torch.nn as nn
import numpy as np

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

class SimpleVMCWavefunction(nn.Module):
    def __init__(self, n_electrons, n_nuclei):
        super().__init__()
        self.n_electrons = n_electrons
        
        # Simple neural network - input: electron-nuclear distances
        input_size = n_electrons * n_nuclei
        self.net = nn.Sequential(
            nn.Linear(input_size, 16),
            nn.Tanh(),
            nn.Linear(16, 8),
            nn.Tanh(),
            nn.Linear(8, 1),
            nn.Tanh()  # Output between -1 and 1
        )
        
        # Very small initialization
        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.001)
                nn.init.zeros_(m.bias)
    
    def forward(self, r, nuclei, charges):
        """
        r: (n_electrons, 3) positions
        nuclei: list of nuclear positions
        charges: list of nuclear charges
        """
        # Base wavefunction: product of exponentials
        psi_base = 1.0
        distances = []
        
        for i in range(self.n_electrons):
            for j, nucleus in enumerate(nuclei):
                dist = torch.norm(r[i] - nucleus)
                distances.append(dist)
                # Use reasonable fixed decay rate
                alpha = charges[j] * 0.8  # Slightly screened
                psi_base *= torch.exp(-alpha * dist)
        
        # Neural network correction based on distances
        distances = torch.stack(distances)
        nn_factor = 1.0 + 0.1 * self.net(distances).squeeze()
        
        return psi_base * nn_factor

def compute_energy(r, wf, nuclei, charges, eps=2e-3):
    """Compute local energy with finite differences"""
    psi = wf(r, nuclei, charges)
    
    # Kinetic energy
    kinetic = 0.0
    for i in range(len(r)):
        for d in range(3):
            r_plus = r.clone()
            r_minus = r.clone()
            r_plus[i, d] += eps
            r_minus[i, d] -= eps
            
            psi_plus = wf(r_plus, nuclei, charges)
            psi_minus = wf(r_minus, nuclei, charges)
            
            kinetic -= 0.5 * (psi_plus - 2*psi + psi_minus) / (eps**2) / psi
    
    # Potential energy
    V = 0.0
    
    # Nuclear repulsion
    for i in range(len(nuclei)):
        for j in range(i+1, len(nuclei)):
            V += charges[i] * charges[j] / torch.norm(nuclei[i] - nuclei[j])
    
    # Electron-nuclear attraction
    for i in range(len(r)):
        for j in range(len(nuclei)):
            V -= charges[j] / torch.norm(r[i] - nuclei[j])
    
    # Electron-electron repulsion
    for i in range(len(r)):
        for j in range(i+1, len(r)):
            V += 1 / torch.norm(r[i] - r[j])
    
    return kinetic + V

def run_vmc(nuclei, charges, n_electrons, n_iter=500):
    """Simple VMC with minimal complexity"""
    # Setup
    nuclei = [n.to(device) for n in nuclei]
    charges = torch.tensor(charges, device=device, dtype=torch.float32)
    
    wf = SimpleVMCWavefunction(n_electrons, len(nuclei)).to(device)
    opt = torch.optim.Adam(wf.parameters(), lr=0.001)
    
    # Initial positions
    r = torch.zeros(n_electrons, 3, device=device)
    for i in range(n_electrons):
        r[i] = nuclei[i % len(nuclei)] + 0.1 * torch.randn(3, device=device)
    
    energies = []
    
    for it in range(n_iter):
        # Collect samples
        E_list = []
        psi_list = []
        
        for _ in range(50):
            # Metropolis move
            r_new = r + 0.3 * torch.randn_like(r)
            
            psi_old = wf(r, nuclei, charges)
            psi_new = wf(r_new, nuclei, charges)
            
            if torch.rand(1, device=device) < (psi_new/psi_old)**2:
                r = r_new
            
            # Energy
            E = compute_energy(r, wf, nuclei, charges)
            psi = wf(r, nuclei, charges)
            
            if torch.isfinite(E) and -10 < E < 10:
                E_list.append(E)
                psi_list.append(psi)
        
        if len(E_list) > 10:
            # Convert to tensors
            E_tensor = torch.stack(E_list)
            psi_tensor = torch.stack(psi_list)
            
            # Weighted average
            weights = psi_tensor**2
            weights = weights / weights.sum()
            E_avg = (weights * E_tensor).sum()
            
            energies.append(E_avg.item())
            
            # Simple loss: just minimize energy
            loss = E_avg
            
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(wf.parameters(), 0.1)
            opt.step()
        
        if it % 100 == 0:
            recent = np.mean(energies[-20:]) if len(energies) > 20 else energies[-1]
            print(f"Iter {it}: E = {E_avg.item():.4f}, recent = {recent:.4f}")
    
    # Return median of last quarter
    final_energies = energies[-len(energies)//4:]
    return np.median(final_energies)

# Run experiments
print("\nMinimal Stable VMC\n")

print("Hydrogen atom:")
E_H = run_vmc([torch.tensor([0.,0.,0.])], [1.], 1, n_iter=400)
print(f"Result: {E_H:.4f} (exact: -0.5000)\n")

print("Helium atom:")
E_He = run_vmc([torch.tensor([0.,0.,0.])], [2.], 2, n_iter=600)
print(f"Result: {E_He:.4f} (exact: -2.9037)\n")

print("H2 molecule:")
E_H2 = run_vmc([torch.tensor([-0.7,0,0]), torch.tensor([0.7,0,0])], [1.,1.], 2, n_iter=600)
print(f"Result: {E_H2:.4f} (exact: -1.1700)")