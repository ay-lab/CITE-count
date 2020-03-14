from setuptools import setup, find_packages
from pathlib import Path


def readme():
    with open('README.md') as f:
        return f.read()

requirements = [
        line.strip() for line in Path('requirements.txt').read_text('utf-8').splitlines()
    ]

setup(
    author="Niu Du",
    author_email='jiyin.dna@gmail.com',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Read and QC totalseq A fastq files and export count ADT matrix, merge with sc RNASeq data and export to scanpy",
    install_requires=requirements,
    license="MIT license",
    long_description=readme(),
    long_description_content_type='text/markdown',
    keywords='cite count',
    name='cite_count',
    packages=find_packages(include=['cite_count', 'cite_count.py','utilities.py']),
    url='https://github.com/ndu-UCSD/cite_count',
    version='0.5.2',
    include_package_data=True,
    zip_safe=False,
)

