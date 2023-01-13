from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='usecases_simairr',
    version='0.1',
    packages=['bias_demo_figure','usecases_simairr'],
    url='',
    license='MIT',
    author='Chakravarthi Kanduri',
    author_email='chakra.kanduri@gmail.com',
    description='',
    include_package_data=True,
    zip_safe=False,
    entry_points={'console_scripts': ['multi_simairr=usecases_simairr.multi_simairr:execute',
                                      'multi_ml=usecases_simairr.multi_ml:execute',
                                      'concat_airr=bias_demo_figure.convert_airr_to_compairr_inputs:execute',
                                      'gen_compairr_input=bias_demo_figure.convert_olga_to_compairr_inputs:execute',
                                      'run_compairr_seqcounts=bias_demo_figure.run_compairr_seqcounts:execute',
                                      'compute_pgen_public=bias_demo_figure.compute_pgen_public_sequences:execute',
                                      'compute_pval=bias_demo_figure.compute_pval:execute',
                                      'concat_pdata=bias_demo_figure.concatenate_pdata:execute']})