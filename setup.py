from setuptools import setup, find_packages

with open('README.md') as readme_file:
    README = readme_file.read()

with open('HISTORY.md') as history_file:
    HISTORY = history_file.read()

setup_args = dict(
    name='PeptideSearch',
    version='0.1',
    description='Search a list of peptides in a fasta file (or proteome) and returns exact matches or hits with one mismatch',
    long_description_content_type="text/markdown",
    long_description=README + '\n\n' + HISTORY,
    license='MIT',
    packages=find_packages(),
    author='Ehdieh Khaledian',
    author_email='khaledianehdieh@gmail.com',
    keywords=['Peptide Search', 'Peptide Match', 'One Mismatch', 'Exact Match', 'python', 'Peptide', 'Protein'],
    url='https://github.com/khaledianehdieh/PeptideSearch',
    download_url='https://pypi.org/project/PeptideSearch/'
)

install_requires = ['pandas', 'numpy', 're', 'itertools', 'more_itertools']

if __name__ == '__main__':
    setup(**setup_args, install_requires=install_requires)
