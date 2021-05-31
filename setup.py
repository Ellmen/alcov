from setuptools import setup

setup(
    name='alcov',
    version='0.1',
    description='Identify frequencies of concerning mutations from aligned reads',
    author='Isaac Ellmen',
    author_email='isaac.ellmen@uwaterloo.ca',
    packages=['alcov'],
    install_requires=[
        'fire',
        'numpy',
        'pandas',
        'scikit-learn>=0.24',
        'matplotlib',
        'seaborn',
        'pysam',
    ],
    entry_points={
        'console_scripts': ['alcov=alcov.command_line:main'],
    }
)
