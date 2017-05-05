import pandas as pd
import settings


def _load_conversion_dict(source, target, id_type):
    fname = settings.RESOURCE_DIR + '/conversion/%s_id_conversion.tsv' % id_type
    conversion_df = pd.read_csv(fname, sep='\t', dtype=str)
    conversion_df = conversion_df[conversion_df[source].notnull() &
                                  conversion_df[target].notnull()]
    return dict(zip(conversion_df[source], conversion_df[target]))

def _convert_list(id_list, conversion_dict, ignore_null):
    if ignore_null:
        return [conversion_dict[str(id)] for id in id_list if id in conversion_dict]
    else:
        return [conversion_dict[str(id)] if id in conversion_dict else None 
                                         for id in id_list]

def convert(id_list, source, target, id_type='gene', ignore_null=False):
    conversion_dict = _load_conversion_dict(source, target, id_type)
    if type(id_list[0]) == type([]):
        converted = []
        for ids in id_list:
            converted.append(_convert_list(ids, conversion_dict, ignore_null))
        return converted
    else:
        return _convert_list(id_list, conversion_dict, ignore_null)
