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


def create_histogram_all(args, mode = "all", samples = None):
    """A function which plots a histogram showing the distribution of a
    chosen element across molecules for all databases.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    element = input(
        """Enter the symbol of the element you would like to create the
        histogram for (if you don't know what elements to expect, see the
        molecular composition first).
        """
        )
    num_datasets = len(args)
    if mode == "compare":
        bar_width = 0.8 / (num_datasets * 2)
        shifts = np.arange(-num_datasets, num_datasets) * bar_width
    else:
        bar_width = 0.8 / num_datasets
        shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width + bar_width/2

    plt.figure()#figsize=(16, 12))

    
    
    for boolean in mode_dict[mode]:
        hatches = ["","","//","","//","","//"]
        for i, arg_set in enumerate(args):
            # if len(args) == 5 and i % 2 != 0:
            #     continue
            if i == 0 or samples == None:
                ecount_arr = get_element_counts(arg_set[1], element, boolean)
            else:
                ecount_arr = get_element_counts(arg_set[1], element, boolean, samples[i-1])
            bins = np.arange(ecount_arr.max() + 1) - 0.5
            shifted_bins = bins + shifts[i]
            label_db=set_label(arg_set[0])
            if boolean == False:# and mode == "compare":
                shifted_bins = bins + shifts[:num_datasets][i]
                label = label_db#f"{arg_set[0]} (all molecules)"
            else:
                label = label_db+"(valid only)"#f"{arg_set[0]} (valid only)"
                if mode == "compare":
                    shifted_bins = bins + shifts[num_datasets:][i]
                else:
                    shifted_bins = bins + shifts[i]

            data = np.array([ecount_arr, shifted_bins], dtype=object)
            np.save(f"/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/data_{element}_{label}", data)


            # n,b,patches =
            plt.hist(
                ecount_arr,
                bins = shifted_bins,
                edgecolor = "black",
                width = bar_width,
                rwidth = bar_width,
                label = label,
                density = True
                )
            # for p in patches:
            #     p.set_hatch(hatches[i])

        plt.xlabel("Count",fontsize=35)
        plt.ylabel("Density",fontsize=35)
        plt.tick_params(axis='both', labelsize=25)

    #plt.xticks(ticks=range(0, 19, 1))  # 19 for h 
    #plt.xlim(left=-0.5)

    # if len(args) == 8:
    #     plt.title(f"Histogram Showing Distribution of {element} in thiols databases")
    # elif len(args) == 2:
    #     plt.title(f"Histogram Showing Distribution of {element} in qm9 databases")
    # elif len(args) == 5:
    #     plt.title(f"Histogram Showing Distribution of {element} in OE62 databases")
    # else:
    #     plt.title(f"Histogram Showing Distribution of {element} in OE62+THz databases")
    plt.margins(x=0)
    plt.legend()
    plt.tight_layout()

    
    plt.show()

    return 

def plot_kde_hist_all(args, mode = "all", samples = None):
    """This function plots the kde on top of a normalised histogram for all
    databases, so that it can be seen that the two distributions are
    consistent.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    element = input(
        """Enter the symbol of the element you would like
        to create the histogram for (if you don't know what elements to
        expect, see the molecular composition first).
        """
        )
    num_datasets = len(args)
    if mode == "compare":
        bar_width = 0.8 / (num_datasets * 2)
        shifts = np.arange(-num_datasets, num_datasets) * bar_width
    else:
        bar_width = 0.8 / num_datasets
        shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width

    plt.figure()#figsize=(16, 12))

    for boolean in mode_dict[mode]:
        for i, arg_set in enumerate(args):
            if len(args) == 5 and i % 2 != 0 and mode != "all":
                continue
            # Plotting histogram
            if i == 0 or samples == None:
                ecount_arr = get_element_counts(arg_set[1], element, boolean)
            else:
                ecount_arr = get_element_counts(arg_set[1], element, boolean, samples[i-1])
            bins = np.arange(ecount_arr.max() + 1) - 0.5
            label_db=set_label(arg_set[0])
            if boolean == False:
                shifted_bins = bins + shifts[:int(num_datasets)][i]
                label = label_db#f"{arg_set[0]} (all molecules)"
            else:
                label = label_db+"(valid)"# f"{arg_set[0]} (valid only)"
                if mode == "compare":
                    shifted_bins = bins + shifts[int(num_datasets):][i]
                else:
                    shifted_bins = bins + shifts[i]
            plt.hist(
                ecount_arr,
                bins = shifted_bins,
                edgecolor = "black",
                width = bar_width,
                rwidth = bar_width,
                label = label,
                density = True
                )

            # Plotting KDE
            density = gaussian_kde(ecount_arr, bw_method=0.3)
            x = np.linspace(0, ecount_arr.max(), 1000)
            plt.plot(x, density(x), label=arg_set[0], alpha=0.7)

    plt.xlabel("Count")
    plt.ylabel("Density")
    # if len(args) == 8:
    #     plt.title(f"Histogram and KDE Showing Distribution of {element} in thiols databases")
    # elif len(args) == 2:
    #     plt.title(f"Histogram and KDE Showing Distribution of {element} in qm9 databases")
    # elif len(args) == 5:
    #     plt.title(f"Histogram and KDE Showing Distribution of {element} in OE62 databases")
    # else:
    #     plt.title(f"Histogram and KDE Showing Distribution of {element} in OE62+THz databases")

    plt.margins(x=0)
    plt.legend()
    plt.tight_layout()
    plt.show()




def get_element_counts(data_file, symbol, use_valid = False, sample = None):
    """A function which gets the count of a chosen element in each molecule."""

    f = open(f"{base_dir}/{data_file}", "r", encoding = "cp1252")
    ecomp_data = json.load(f)
    if sample is not None:
        summary = ecomp_data[-1]
        ecomp_data = [ecomp for i, ecomp in enumerate(ecomp_data) if i in sample]
        ecomp_data.append(summary)
    f.close()

    

    molecule_count = len(ecomp_data) - 1
    ecount_arr = np.zeros((molecule_count,), dtype=int)

    for molecule_num in range(0, molecule_count):
        try:
            ecount_arr[molecule_num] = ecomp_data[molecule_num][symbol]
        except KeyError:
            pass

    return ecount_arr


def plot_atoms_data_all(args, mode = "all", samples = None, split = False):
    """Plots distributions for number of atoms and molecular weights in all
    databases.
    """
    homo_lumo=False
    split_dbs = False
    if split == True and len(args) == 3:
        args[0][0] = "OE62"
        args.insert(1, ("THz", args[0][1]))
        args[2][0] = "Generated\n(transform)"
        args[3][0] = "Generated\n(no transform)"
        split_dbs = True
    num_datasets = len(args)
    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    plt.figure(figsize=(18, 8))
    no_of_atoms_list = []
    weights_list = []
    hl_list = []
    x_lims = [np.inf, 0]
    y_lims = [np.inf, 0]

    for boolean in mode_dict[mode]:
        for i, arg_set in enumerate(args):
            if (len(args) == 5 and i % 2 != 0 and mode != "all" and split_dbs == False) or (split_dbs == True and i == 1):
                split = False
                continue
            if ((i == 0 or samples == None) and split == False):
                if homo_lumo:
                    no_of_atoms, weights,hl = get_atoms_data(arg_set[1], boolean, split = split,homo_lumo=homo_lumo)
                    hl_list.append(hl)                
                else:
                    no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, split = split)
                
                no_of_atoms_list.append(no_of_atoms)
                weights_list.append(weights)
                
            elif split == True and i == 0:
                no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, split = split)
                for i in range(0, 2):
                    if samples != None:
                        no_of_atoms[i] = [atoms for j, atoms in enumerate(no_of_atoms[i]) if j in samples[i]]
                        weights[i] = [weight for j, weight in enumerate(weights[i]) if j in samples[i]]
                    no_of_atoms_list.append(no_of_atoms[i])
                    weights_list.append(weights[i])

                    if max(no_of_atoms[i]) > x_lims[1]:
                        x_lims[1] = max(no_of_atoms[i])
                    if max(weights[i]) > y_lims[1]:
                        y_lims[1] = max(weights[i])
                    if min(no_of_atoms[i]) < x_lims[0]:
                        x_lims[0] = min(no_of_atoms[i])
                    if min(weights[i]) < y_lims[0]:
                        y_lims[0] = min(weights[i])
            else:
                if split_dbs == True:
                    no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, sample = samples[i], split = split)
                elif homo_lumo:
                    no_of_atoms, weights, hl = get_atoms_data(arg_set[1], boolean, sample = samples[i - 1], split = split,homo_lumo=homo_lumo)
                    hl_list.append(hl)
                    print("HOMO LUMO ", np.shape(weights_list), np.shape(hl_list))
                else:
                    no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, sample = samples[i - 1], split = split)

                no_of_atoms_list.append(no_of_atoms)
                weights_list.append(weights)
                
                
                
                if max(no_of_atoms) > x_lims[1]:
                    x_lims[1] = max(no_of_atoms)
                if max(weights) > y_lims[1]:
                    y_lims[1] = max(weights)
                if min(no_of_atoms) < x_lims[0]:
                    x_lims[0] = min(no_of_atoms)
                if min(weights) < y_lims[0]:
                    y_lims[0] = min(weights)

        data = np.array([no_of_atoms_list],dtype=object)
        np.save(f"/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/no_of_atoms_unfiltered", data)

        for i, arg_set in enumerate(args):
            if (len(args) == 5 and i % 2 != 0 and mode != "all"):
                continue
            density = gaussian_kde(no_of_atoms_list[i], bw_method=0.3)
            x = np.linspace(0, max(no_of_atoms_list[i]), 1000)
            label_db=set_label(arg_set[0])
            if boolean == True:
                label = label_db#f"{arg_set[0]} (valid only)"
            else:
                label = label_db #f"{arg_set[0]} (all molecules)"
            plt.plot(x, density(x), label=label, alpha=0.7)
            print("mean value: ", x[np.argmax(density(x))], np.mean(no_of_atoms_list[i]))

    plt.xlabel("Count", fontsize = 40)
    plt.ylabel("Density", fontsize = 40)
    # if len(args) == 8:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in thiols databases")
    # elif len(args) == 2:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in qm9 databases")
    # elif len(args) == 5:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in OE62 databases")
    # else:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in OE62+THz databases")
    
    plt.margins(x=0)
    plt.xlim(0,70)#loc=[1.05,0])
    plt.legend(fontsize=30)#(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.tight_layout()
    plt.show()

    if len(args) == 5 and mode != "all":
        args = [args[i] for i in range(0, len(args)) if i % 2 == 0]
    if mode == "compare":
        args = args * 2

    plt.figure(figsize=(18, 8))

    for i, arg_set in enumerate(args):
        density = gaussian_kde(weights_list[i], bw_method=0.3)
        x = np.linspace(0, max(weights_list[i]), 1000)
        label_db=set_label(arg_set[0])
        if boolean == True:
            label = label_db #f"{arg_set[0]} (valid only)"
        else:
            label = label_db #f"{arg_set[0]} (all molecules)"
        plt.plot(x, density(x), label=label, alpha=0.7)

    plt.xlabel("Weight (g/mol)", fontsize = 40)
    plt.ylabel("Density", fontsize = 40)
    # if num_datasets == 8:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in thiols databases")
    # elif num_datasets == 2:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in qm9 databases")
    # elif num_datasets == 5:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62 databases")
    # else:
        # plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62+THz databases")
    plt.xlim(70,1000)#170)#
    plt.margins(x=0)
    plt.legend(fontsize=30)#(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.tight_layout()
    plt.savefig(f"/root/MChem_DGMs/analysis/Plots/DRUGS/weights.pdf",format='pdf',dpi=300)
    plt.show()

   # plot HOMO-LUMO gap prediction
    if homo_lumo:
        for i, arg_set in enumerate(args):
            density = gaussian_kde(hl_list[i], bw_method=0.3)
            x = np.linspace(0, max(hl_list[i]), 1000)
            label_db=set_label(arg_set[0])
            if boolean == True:
                label = label_db #f"{arg_set[0]} (valid only)"
            else:
                label = label_db #f"{arg_set[0]} (all molecules)"
            plt.plot(x, density(x), label=label, alpha=0.7)

        plt.xlabel("HOMO-LUMO gap /eV")
        plt.ylabel("Density")
        # if num_datasets == 8:
        #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in thiols databases")
        # elif num_datasets == 2:
        #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in qm9 databases")
        # elif num_datasets == 5:
        #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62 databases")
        # else:
            # plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62+THz databases")
        plt.xlim(0,17.5)#170)#
        plt.margins(x=0)
        plt.legend()
        plt.tight_layout()
        plt.show()


    # Creating a 2D density plot with accompanying histograms on respective
    # axes for each database.
    for i, arg_set in enumerate(args):
        arg_set = list(arg_set)
        label_db=set_label(arg_set[0])
        if (mode == "compare" and i >= len(args)/2) or mode == "valid":
            arg_set[0] = label_db #f"{arg_set[0]} (valid only)"
        else:
            arg_set[0] = label_db #f"{arg_set[0]} (all molecules)"
        # plot_2D_density_atom_data(arg_set[0], no_of_atoms_list[i], weights_list[i], x_lims, y_lims)

def get_atoms_data(data_file, use_valid, sample = None, split = False, homo_lumo=False):
    """A function which gets the number of atoms in, and molecular weight 
    of, each molecule in a database.
    """

    if use_valid == False:
        f = open(f"{base_dir}/{data_file}", "r", encoding = "cp1252")
        ecomp_data = json.load(f)
        if sample is not None:
            summary = ecomp_data[-1]
            ecomp_data = [ecomp for i, ecomp in enumerate(ecomp_data) if i in sample]
            ecomp_data.append(summary)
        f.close()
        # loading in homo-lumo gap predictions for oe62
        if homo_lumo:
            if "of_" in data_file:
                hl=np.load("/home/zkoczor//Work2/cgschnet-thiols/gschnet/OE62/homo_lumo/original_training/homo_lumo.npz",allow_pickle=True)["HL"][()]
    
            if "o2_" in data_file:
                hl=np.load("/home/zkoczor//Work2/cgschnet-thiols/gschnet/OE62/homo_lumo/OE62_2_filtered_generated/homo_lumo.npz",allow_pickle=True)["HL"][()]
    # print(min(hl),max(hl)) 
    # hl_train=np.load("/home/zkoczor//Work2/cgschnet-thiols/gschnet/OE62/homo_lumo/original_training/homo_lumo.npz",allow_pickle=True)["HL"][()]

            if sample is not None:
                hl = [ecomp for i, ecomp in enumerate(hl) if i in sample]
    else:
        ecomp_data = get_valid(data_file, sample)

    if split == True:
        no_of_atoms = [[], []]
        weights = [[], []]
        for molecule in ecomp_data[:-1]:
            if "Au" in molecule:
                no_of_atoms[1].append(sum(list(molecule.values())))
                molecular_weight = 0
                for symbol, ecount in molecule.items():
                    molecular_weight += atomic_masses[atomic_numbers[symbol]] * ecount
                weights[1].append(molecular_weight)
            else:
                no_of_atoms[0].append(sum(list(molecule.values())))
                molecular_weight = 0
                for symbol, ecount in molecule.items():
                    molecular_weight += atomic_masses[atomic_numbers[symbol]] * ecount
                weights[0].append(molecular_weight)

    else:    
        no_of_atoms = [sum(list(molecule.values())) for molecule in ecomp_data[:-1]]
        weights = []
        for molecule in ecomp_data[:-1]:
            molecular_weight = 0
            for symbol, ecount in molecule.items():
                molecular_weight += atomic_masses[atomic_numbers[symbol]] * ecount
            weights.append(molecular_weight)

    if homo_lumo:
        return no_of_atoms, weights, hl
    return no_of_atoms, weights