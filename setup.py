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
    raise RuntimeError("Unable to find version string.")


with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

# now that pyproject.toml is used, we should only need a minimal setup.py
# to ensure compatibility with older versions of pip
# see https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html
# the full setup.py is commented out below, will be removed in future.
setup()

# setup(
#     name='ampere',
#     version=get_version("ampere/__init__.py"),
#     description='',
#     long_description=long_description,
#     long_description_content_type="text/markdown",
#     url='https://github.com/icsm/ampere',
#     project_urls={
#         "Bug Tracker": "https://github.com/icsm/ampere/issues",
#     },
#     author=('Peter Scicluna, Francisca Kemper, Sundar Srinivasan, Jonathan'
#             ' Marshall, Oscar Morata, Alfonso Trejo, Sascha Zeegers, Lapo'
#             ' Fanciullo, Thavisha Dharmawardena'),
#     author_email='peter.scicluna@eso.org',
#     license='GPL',
#     packages=['ampere',
#               'ampere.infer',
#               'ampere.utils',
#               'ampere.data',
#               'ampere.models'],
#     package_data={'ampere': ['ampere_allfilters.hd5', ], },
#     install_requires=['numpy',
#                       'astropy>=4.0.0',
#                       'scipy',
#                       'matplotlib',
#                       'spectres',
#                       'tqdm',
#                       'pyphot',
#                       'emcee',
#                       'dynesty',
#                       'corner',
#                       ],
#     python_requires=">=3.7",
#     extras_require={
#         "sbi": ["torch", "sbi"],
#         "zeus": ["zeus-mcmc"],
#         "arviz": ["arviz"],
#         "dev": ['sphinx',  # anything for development, tests or docs
#                 'nbsphinx',
#                 'nbconvert'],
#         },
#     classifiers=[
#         'Development Status :: 3 - Alpha',
#         'Intended Audience :: Science/Research',
#         'License :: OSI Approved :: GPL License',
#         'Operating System :: POSIX :: Linux',
#         'Programming Language :: Python :: 3',
#     ],
# )
