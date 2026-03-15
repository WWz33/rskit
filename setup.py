from setuptools import setup, find_packages

setup(
    name="rnaseq-tools",
    version="0.1.0",
    description="Python interface for RNA-seq analysis tools (STAR + Salmon)",
    packages=find_packages(),
    python_requires=">=3.8",
    entry_points={"console_scripts": ["rnaseq-tools=rnaseq_tools.cli:main"]},
)
