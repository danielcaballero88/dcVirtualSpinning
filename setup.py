from setuptools import setup
import os 

here = os.path.dirname(".")

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    # Basic setup configuration
    name='dcVirtualSpinning',
    version='0.1dev',
    description='Virtual (Electro)Spinning Microstructure Generator',
    long_description=long_description,  # same as README.md
    long_description_content_type='text/markdown',
    url='https://github.com/khawabonga/dcVirtualSpinning',
    author='khawabonga',
    author_email='khawabonga@gmail.com',
    #
    # Classifiers help users find your project by categorizing it.
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 3',
    ],
    keywords='multiscale simulation',
    packages=['VirtualSpinning'], # se podria hacer = find_packages()
    python_requires='>3.6',
    license='MIT',
)
