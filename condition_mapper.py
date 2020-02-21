def map_inputs(user, all_data, data_name, id_list):
    '''
    This maps the array of group identifiers from the front-end
    and places them in an array

    The front-end sends back a data structure like:
    {"0": {"0":"group1", "1":"group2"}, "1": {"0":"groupA", "1":"groupB"}}
    so that in row 0 we have a comparison of group2 vs group1, and so on.
  
    id_list has the names of the WDL inputs.
    '''
    unmapped_data = all_data[data_name]
    
    col_num_to_id_mapping = dict(zip(range(len(id_list)), id_list))

    d = {}
    for i in id_list:
        d[i] = []

    for row_num, row_dict in unmapped_data.items():
        for col_id, val in row_dict.items():
            col_num = int(col_id) # e.g. "1" --> 1
            target_var = col_num_to_id_mapping[col_num] # e.g. exchange 1 for "Base"
            d[target_var].append(val)
    return d
