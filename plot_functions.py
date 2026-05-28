## for making smiles.json

def mol2smiles(mol):
    try:
        Chem.SanitizeMol(mol)
    except ValueError:
        return None
    return Chem.MolToSmiles(mol)

## for making ecomps.json

def generate_data(database, data_file):
    """Defining a function to generate elemental composition analysis data
    for a chosen database.
    """
    # print(sqlite3.version) 
    # print(sqlite3.sqlite_version)
    # Connecting to one of the databases.

    db = ase.db.connect(f"{base_dir}/{database}")

    # Initialising an empty dictionary to store the number of times each
    # element occurs in total.
    element_counts_total = {}

    #Initialising an empty list to store all data objects.
    all_element_data = []

    # Iterating through the database to get the number of atoms of each
    # element for each molecule.
    for row in db.select():
        atoms = row.toatoms()
        elements = atoms.get_chemical_symbols()
        # Initialising an empty dictionary to store the number of times each
        # element occurs in each molecule (resets every iteration)
        element_counts = {}

        # Another loop which adds the count of each element to the dictionaries.
        for i in set(elements):
            element_count = elements.count(i)
            element_counts[i] = element_count

            if i not in element_counts_total:
                element_counts_total[i] = element_count

            else:
                element_counts_total[i] += element_count

        # Adds the dictionaries containing elemental composition data for each
        # molecule to the empty list.
        all_element_data.append(element_counts)

    # Adds the dictionary containing elemental composition data for the
    # whole database to the empty list, then dumps the list to a json file.
    all_element_data.append(element_counts_total)
    f = open(f"{base_dir}/{data_file}", "w", encoding = "cp1252")
    json.dump(all_element_data, f)
    f.close()



