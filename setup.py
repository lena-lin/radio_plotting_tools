from setuptools import setup, find_packages

setup(
    name='radio_plotting_tools',
    author='Lena Linhoff, Kevin Schmidt',
    author_email='lena.linhoff@tu-dortmund.de',
    version='0.0.1',
    packages=find_packages(),
    install_requires=[
        'astropy',
        'numpy',
        'matplotlib',
        'click',
    ],
    entry_points={
        'console_scripts': [
            'radio_plotting_clean_map = radio_plotting_tools.scripts.plot_clean_map:main',
            'make_gif = radio_plotting_tools.scripts.make_gif:main',
        ]
    }
)
