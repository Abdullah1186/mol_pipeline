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

def get_average_composition(data_file, sample = None):
    """A function which takes the dictionaries containing the number of
    atoms of each element in each molecule, finds the percentage elemental
    composition of each and divides it by the number of molecules. The 
    result is displayed to the user.
    """

    
    f = open(f"{base_dir}/{data_file}", "r", encoding = "cp1252")
    ecomp_data = json.load(f)
    if sample is not None:
        summary = ecomp_data[-1]
        ecomp_data = [ecomp for i, ecomp in enumerate(ecomp_data) if i in sample]
        ecomp_data.append(summary)
    f.close()
    
    avg_mol_comp = {}
    

    for molecule in ecomp_data:
        
        num_atoms = sum(list(molecule.values()))
        for symbol, count in molecule.items():
            try:
                avg_mol_comp[symbol] += count/num_atoms
            except KeyError:
                avg_mol_comp[symbol] = count/num_atoms
    
    
    for element in avg_mol_comp:
        avg_mol_comp[element] /= len(ecomp_data) - 1
    
    return avg_mol_comp

def get_average_composition_num(data_file,  sample = None, ):
    """A function which converts the symbol keys in the average molecular
    composition dictionary to atomic numbers.
    atomic numbers come from ase 
    """

    avg_mol_comp = get_average_composition(data_file, sample)
    
    avg_mol_comp_num = {}
    for element, count in avg_mol_comp.items():
        atomic_num = atomic_numbers[element]
        avg_mol_comp_num[atomic_num] = count

    return avg_mol_comp_num

def create_bar_all(args, mode = "all", samples = None, split = False):
    """A function which generates a barplot showing average atom percentage
    per molecule in all databases. This is done in terms of atomic number
    so it is easier to compare average size of molecules between
    databases. the mode is always all never compare and split is always false 
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    
    num_datasets = len(args)
    all_elements = []
    avg_mol_comp_list = []
    atomic_nums_list = []

    # The purpose of this loop is primarily to get a list of all of the
    # elements across all of the databases. Lists of other data like
    # the dictionary of elements to average molecular composition and lists
    # of atomic numbers present in each database are obtained along the way.
    for boolean in mode_dict[mode]:
        for i, arg_set in enumerate(args):
            if i == 0 or samples == None:
                avg_mol_comp_num = get_average_composition_num(arg_set[1], boolean, split = split)
            else:
                avg_mol_comp_num = get_average_composition_num(arg_set[1], boolean, samples[i-1], split = split)

            
            atomic_nums_list.append(list(map(str, sorted(list(avg_mol_comp_num.keys())))))
            avg_mol_comp_list.append(avg_mol_comp_num)
            if set(atomic_nums_list[-1]) != set(all_elements):
                all_elements += list(set(atomic_nums_list[-1]).difference(set(all_elements)))

        # The list of all atomic numbers is sorted as though they were integer values.
        all_elements = sorted(all_elements, key = int)

        # If there are any elements which exist in one of the databases but
        # not another, a key is created and set to 0 where it is absent.
        an_index = 0
        for atomic_nums in atomic_nums_list:
            for element in all_elements:
                if element not in atomic_nums:
                    avg_mol_comp_list[an_index][int(element)] = 0
            an_index += 1

        # x and y values obtained
        average_counts = list(map(lambda x:list(dict(sorted(x.items())).values()), avg_mol_comp_list))
        num_av_counts = np.array(average_counts)
        x = np.arange(len(all_elements))   
        # np.save("my_array",num_av_counts)
        # np.save("x_arr",x)


        # Adjusting the positions of the different plots
        # so they are all visible and changing their size so they all fit.
        plt.figure(figsize=(9,5.5))
        # N=9
        # plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.RdYlGn(np.linspace(0,1,N)))

        bar_width = 0.1#0.3
        if mode == "compare":
            shifts = np.arange(-num_datasets, num_datasets) * bar_width + bar_width/2
        else:
            shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width + bar_width/2
        if not(mode == "compare" and boolean == True):
            shifted_xs = [x + shifts[:num_datasets][i] for i in range(num_datasets)]
        else:
            shifted_xs = [x + shifts[num_datasets:][i] for i in range(num_datasets)]

        # Finally, the plot is made.
        for i in range(num_datasets):
            j = i
            label_db=set_label(args[i][0])
            if boolean == True and mode == "compare":
                j += num_datasets
                label = label_db #f"{args[i][0]} (valid only)"
            else:
                label = label_db #f"{args[i][0]} (all molecules)"
            plt.bar(
                shifted_xs[i],
                average_counts[j],
                bar_width,
                label = label,
                edgecolor = "black"
                )

    plt.xlabel('Atomic Number',fontsize=35)
    plt.ylabel('Average Proportion per Molecule',fontsize=30)
    

    plt.xticks(np.arange(len(all_elements)), labels = all_elements, fontsize=35)
    plt.yticks(fontsize=35)
    plt.margins(x=0)
    plt.legend(fontsize=35)
    plt.tight_layout()
    plt.show()
    return x, num_av_counts

