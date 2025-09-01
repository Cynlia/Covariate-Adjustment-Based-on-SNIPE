
# Covariate Adjustment Cannot Hurt: Treatment Effect Estimation under Network Interference with Low-Order Outcome Interactions

[![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b.svg)](https://arxiv.org/abs/XXXX.XXXXX)

This repository contains the code for our paper:

**Covariate Adjustment Cannot Hurt: Treatment Effect Estimation under Interference with Low-Order Outcome Interactions**  
Xinyi Wang, Shuangning Li  
<!--Published at [Conference/Journal Name], [Year]-->  
[Link to arXiv or publication]

---

## 📦 Repository Structure

```
.
├── README.md                # Project overview
├── LICENSE                  # License file
├── data/                    # Data scripts or instructions
│   └── README.md
├── src/                     # Core method implementations
│   ├── Codes_for_Experiments/
│       ├── experiments_master_graph_aware.py
│       ├── nci_linear_setup.py
│       └── nci_polynomial_setup.py
│   ├── Codes_for_Plots/
│       ├── master_plots_graph_aware_MSE.py
│       └── master_plots_graph_aware.py
│   └── OutFiles/.           # Output files (figures, csvs)
│       ├── graph_aware/.    # Outputs in paper
│       └── new/.                        
├── notebooks/               # Jupyter notebooks for demo
│   └── demo.ipynb
└── paper/                   # Copy of the paper
    └── your_paper.pdf
```

---

## 🚀 Clone the repo

```bash
git clone https://github.com/Cynlia/Covariate-Adjustment-Based-on-SNIPE.git
cd Covariate-Adjustment-Based-on-SNIPE
```
---

## ▶️ Running the Code

To run the main experiment:

```bash
python src/Code_for_Experiments/experiments_master_graph_aware.py
```

---

## 🙏 Code Attribution

This project is developed based on the following open-source repository:

- [mayscortez/low-order-unitRD](https://github.com/mayscortez/low-order-unitRD):
  We adapted some components of their implementation, especially the experimental design and data generation. Several files in `src/` are modified versions of their code to fit the objectives of our study.

We sincerely thank the original authors for making their code publicly available. Please refer to the respective repositories for licensing terms.

---

## 📄 License

This code is released under the [MIT License](LICENSE), except for components adapted from third-party sources, which retain their original licenses.

---

## 📝 Citation

If you use this repository, please cite our paper:

```bibtex
@inproceedings{yourbibtexentry,
  title={Your Paper Title},
  author={Your Name and Collaborator},
  year={2025},
  journal={Conference/Journal Name}
}
```

---

## 📬 Contact

For questions or collaborations, please contact:  
Xinyi Wang – wang.xinyi@berkeley.edu  
Shuangning Li - shuangning.li@chicagobooth.edu

---
