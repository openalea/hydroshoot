from os.path import isfile


def verify_inputs(path_project: str, **kwargs):
    if 'psi_soil' not in kwargs:
        assert (isfile(f'{path_project}psi_soil.input')), "The 'psi_soil.input' file is missing."
