from setuptools import setup, find_packages

setup(
    name="aria-omics",
    version="4.3.13",
    description="Agentic Research Intelligence for -omics Analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.10",
    packages=find_packages(),
    install_requires=[
        "anthropic>=0.25.0",
        "litellm>=1.30.0",
        "rich>=13.7.0",
        "tiktoken>=0.6.0",
        "pyyaml>=6.0",
        "requests>=2.31.0",
        "tqdm>=4.66.0",
        "click>=8.1.0",
        "numpy>=1.24.0",
        "pandas>=2.0.0",
        "h5py>=3.9.0",
        "anndata>=0.10.0",
    ],
    extras_require={
        "analysis": [
            "scanpy[leiden]>=1.9.6",
            "squidpy>=1.3.0",
            "scipy>=1.11.0",
            "scikit-learn>=1.3.0",
            "matplotlib>=3.7.0",
            "seaborn>=0.12.0",
        ],
        "dev": [
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",
            "black>=23.0.0",
            "ruff>=0.1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "aria=aria.tui:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
    ],
)
