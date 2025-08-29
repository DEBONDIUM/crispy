# Crispy
# 🧩 Crispy – Fragments Analysis in Solid Mechanics

**Crispy** is a lightweight Python library designed to **identify, count, and analyze fragments** generated in **particule-based simulations**. Whether you're working on dynamic fracture, debris tracking, or material fragmentation, `crispy` provides simple and robust tools to extract meaningful metrics from your numerical data.

> 💥 From breakage to data — fast, clear, and **crispy**.

---

## 📦 Features

- 🔍 Detect and label fragments from discrete simulation outputs
- 📊 Compute fragment size distributions and statistics
- 🧱 Support for 2D and 3D domains (nodal position and connectivity)
- 🧮 Works with numpy arrays, mesh files, or custom simulation outputs
- 🧰 Minimal dependencies, easy to integrate into existing workflows

---

## 🚀 Getting Started

### Installation

```bash
pip install git+https://github.com/DEBONDIUM/crispy.git
```

You can also clone and install manually:

```bash
git clone https://github.com/DEBONDIUM/crispy.git
cd crispy
pip install .
```


## 🧪 Example Usage

```python
import crispy as cp

detector = cp.FragmentDetector("path/to/mesh/files")
detector.build_fragments()

for i in range(detector.iterations_nb):
    detector.plot2D(iteration=i, save=True)
```

More examples are available in the [`examples/`](examples/) folder and in the [documentation](https://debondium.github.io/crispy).

---

## 📘 Documentation

The full documentation (installation, API reference, usage examples) is available at:  
🔗 https://debondium.github.io/crispy

---

## 🔬 Applications

- Fracture mechanics  
- Impact simulations  
- Granular material studies  
- Discrete Element Method (DEM) post-processing  
- Crack propagation analysis  

---

## 🛠 Dependencies

- [Numpy](https://numpy.org)  
- [Scipy](https://scipy.org)  
- [Open3d](http://www.open3d.org)  
- [Matplotlib](https://matplotlib.org)

---

## 🧑‍💻 Contributing

Pull requests are welcome!  
Feel free to fork the repository, propose new features, or report bugs via GitHub Issues.

---

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

---

## 👨‍🏫 Acknowledgments

This library is developed for researchers and engineers working in **solid mechanics**, **computational physics**, and **fracture modeling**.  
If you use `crispy` in your work, please consider citing it or linking to the repository.

> _“In the breaking of things lies the story of how they were made.”_
