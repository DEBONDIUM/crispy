# Crispy
# 🧩 Crispy – Fragmentation Analysis in Solid Mechanics

**Crispy** is a lightweight Python library designed to **identify, count, and analyze fragments** generated in **discrete solid mechanics simulations**. Whether you're working on dynamic fracture, debris tracking, or material rupture, `crispy` provides simple and robust tools to extract meaningful metrics from your numerical data.

> 💥 From breakage to data — fast, clear, and **crispy**.

---

## 📦 Features

- 🔍 Detect and label fragments from discrete simulation outputs
- 📊 Compute fragment size distributions and statistics
- 🧱 Support for 2D and 3D domains (structured or unstructured)
- 🧮 Works with numpy arrays, mesh files, or custom simulation outputs
- 🧰 Minimal dependencies, easy to integrate into existing workflows

---

## 🚀 Getting Started

### Installation

```bash
pip install git+https://github.com/DEBONDIUM/crispy.git
```

> Note: You can also clone and install manually:

```bash
git clone https://github.com/DEBONDIUM/crispy.git
cd crispy
pip install .
```

### Example Usage

```python
import crispy

# Load simulation output (e.g., binary mask of fractured domain)
binary_field = load_my_simulation_field()

# Identify fragments
labels = crispy.identify_fragments(binary_field)

# Compute stats
stats = crispy.compute_statistics(labels)

# Visualize
crispy.plot_fragments(labels)
```

---

## 📘 Documentation

Full documentation with examples and API reference is available in the `docs/` folder or at:

📎 [https://your-github-username.github.io/crispy](https://your-github-username.github.io/crispy)

---

## 🔬 Applications

- Fracture mechanics
- Impact simulations
- Granular material studies
- Discrete Element Method (DEM) post-processing
- Crack propagation analysis

---

## 🛠 Dependencies

- `numpy`
- `scipy`
- `matplotlib` (optional, for visualization)

---

## 🧑‍💻 Contributing

Pull requests are welcome! Feel free to fork the repository, propose new features, or report bugs via GitHub Issues.

---

## 📄 License

MIT License — see [LICENSE](LICENSE) file for details.

---

## 👨‍🏫 Acknowledgments

This library is developed for researchers and engineers working in **solid mechanics**, **computational physics**, and **fracture modeling**. If you use `crispy` in your work, please consider citing it or linking to the repository.

---

> _“In the breaking of things lies the story of how they were made.”_
