import os
import sys
from setuptools import setup, Extension, find_packages

if sys.platform == 'win32':
    extra_compile_args = ['/std:c17', '/utf-8']
else:
    extra_compile_args = ['-std=c11', '-finput-charset=UTF-8', '-fexec-charset=UTF-8']

module = Extension(
    'oiiaio_web_kepler',
    sources=['oiiaio_web_kepler.c'],
    libraries=[],
    py_limited_api=True,
    extra_compile_args=extra_compile_args
)

setup(
    name='oiiaio-web-kepler',
    version='1.3.1',
    description='轨道计算器',
    long_description=open('README.md', encoding='utf-8-sig').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    ext_modules=[module],
    install_requires=[],
    python_requires='>=3.8',
    url='https://github.com/oiiaio-web/oiiaio-web-kepler',
    project_urls={
        'Bug Reports': 'https://github.com/oiiaio-web/oiiaio-web-kepler/issues',
        'Source': 'https://github.com/oiiaio-web/oiiaio-web-kepler',
    },
    license='Apache 2.0',
    author='oiiaio-web',
    author_email='3926580760@qq.com',
    classifiers=[
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    options={
        'bdist_wheel': {
            'py_limited_api': 'cp312'
        }
    }
)