import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='facetsAPI',
    version='0.5.0',
    author='Adam Price',
    author_email='pricea2@mskcc.org',
    description='Toolkit for interacting and manipulating Impact FACETS data at MSK.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/mskcc/facetsAPI',
    project_urls = {
        "Bug Tracker": "https://github.com/mskcc/facetsAPI/issues"
    },
    license='MIT',
    packages=['facetsAPI'],
    install_requires=['pandas'],
)
