from . import settings
import os

def _get_all_files_recursively(path):
    all_files = []
    for subpath in os.listdir(path):
        subpath = os.path.join(path, subpath)
        if os.path.isdir(subpath):
            all_files += _get_all_files_recursively(subpath)
        else:
            all_files.append(subpath)
    return all_files

def read_gmt(fl):
    gene_sets = {}
    for line in fl:
        fields = line.strip().split('\t')
        gene_sets[fields[0]] = fields[2:]
    return gene_sets

def read_gene_sets(path):
    fl = open(path)
    extension = fl.name.split('.')[-1].lower()
    if extension == 'gmt':
        return read_gmt(fl)
    else:
        raise ValueError('File %s has unsupported extension %s' % (path, extension))

def read_resource_list(resource_list):
    gene_sets = {}
    for resource in resource_list:
        path = settings.RESOURCE_DIR + resource
        if os.path.isdir(path):
            resource_files = _get_all_files_recursively(path)
        else:
            resource_files = [path]
        for file in resource_files:
            _curr_gene_sets = read_gene_sets(file)
            gene_sets.update(_curr_gene_sets)
    return gene_sets
