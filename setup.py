from setuptools import setup

setup(
    name='cov_breakdown',
    version='0.1',
    description='Identify frequencies of concerning mutations from aligned reads',
    author='Isaac Ellmen',
    author_email='isaac.ellmen@uwaterloo.ca',
    packages=['cov_breakdown'],
    install_requires=[
        'fire',
    ],
    entry_points={
        'console_scripts': ['cov_breakdown=cov_breakdown.command_line:main'],
    }
)
