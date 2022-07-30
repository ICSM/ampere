from setuptools import setup
import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='ampere',
    version=get_version("ampere/__init__.py"),    
    description='',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/icsm/ampere',
    project_urls={
        "Bug Tracker": "https://github.com/icsm/ampere/issues",
    },
    author='Peter Scicluna, Francisca Kemper, Sundar Srinivasan, Jonathan Marshall, Oscar Morata, Alfonso Trejo, Sascha Zeegers, Lapo Fanciullo, Thavisha Dharmawardena',
    author_email='peter.scicluna@eso.org',
    license='GPL',
    packages=['ampere',
              'ampere.infer',
              'ampere.utils',
              'ampere.data',
              'ampere.models'],
    install_requires=['numpy',
                      'astropy>=4.0.0',
                      'scipy',
                      'matplotlib',
                      'spectres',
                      'tqdm',
                      'pyphot',
                      'emcee',
                      'dynesty',
                      'corner',
                      ],
    python_requires=">=3.7",
    extras_require{
        "sbi":["torch", "sbi"],
        "zeus":["zeus-mcmc"],
        "arviz":["arviz"],
        }

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GPL License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
    ],
)
