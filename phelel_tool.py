""" A module full of classes to facilitate Transport Calculations """

import subprocess
import json
import time
from pathlib import Path
from sys import argv

class FileValidation():
    """Essential file and directory validation"""

    def __init__(self, file_path: str):
        """Supply a path to the file validator when initializing."""
        
        self.str_path: str = file_path
        """The path the class will validate."""
        
        self.path_obj: Path = Path(file_path)
        """The Path object version of the supplied path."""
    
    def __eq__(self, op):
        return self.path_obj == op.path_obj

    def __ne__(self, op):
        return self.path_obj != op.path_obj

    def change_path(self, file_path: str):
        """Change the file path the validator is using."""
        self.str_path: str = file_path
        self.path_obj: Path = Path(file_path)

    def __is_valid_path(self) -> bool:
        """Ensure the given path exists on the system."""
        # expand ~ if it exists
        if '~' in self.str_path:
            # If a home directory can’t be resolved, RuntimeError is raised.
            self.path_obj = self.path_obj.expanduser()
        
        # can throw OS error
        self.path_obj = self.path_obj.resolve()

        return self.path_obj.exists()
    
    def valid_file(self) -> Path:
        """Ensure the given path is a valid file on the system."""
        if not self.__is_valid_path() or not self.path_obj.is_file():
            return None
        return self.path_obj
    
    def valid_directory(self) -> Path:
        """Ensure the given path is a valid file on the system."""
        if not self.__is_valid_path() or not self.path_obj.is_dir():
            return None
        return self.path_obj

class VaspInit():
    def __init__(self):
        self.file_validator: FileValidation
        self.poscar_file: Path = None
        self.potpawpbte_dir: Path = None
        self.incar_file: Path = None
        self.vasp_bin_file: Path = None
        self.hdf5_lib_dir: Path = None

        self.poscar_info: dict = None
        """Contains data about the given POSCAR file"""
    
    def copy_init(self, vasp_init):
        """Copy constructor"""
        self.file_validator = vasp_init.file_validator
        self.poscar_file: vasp_init.poscar_file
        self.potpawpbte_dir: vasp_init.potpawpbte_dir
        self.incar_file: vasp_init.incar_file
        self.vasp_bin_file: vasp_init.vasp_bin_file
        self.hdf5_lib_dir: vasp_init.hdf5_lib_dir

        self.poscar_info = vasp_init.poscar_info
    
    
    def __get_input_files(self) -> dict:
        """Get a dictionary of the input files"""
        init_dict = {}

        # check if init.json exists
        init_json = Path(argv[0] + "/../default-input-files/init.json")

        if init_json.exists():
            with init_json.open('r') as f:
                init_dict = json.load(f)
        else:
            raise SystemError("No JSON file found, please perform initialization of the class.")
        
        return init_dict
    
    def __get_poscar_info(self):
        """Get the poscar information into the poscar_info dictionary"""
        
        # TODO: Throw exception
        if not self.poscar_file:
            raise SystemError("Please load in the POSCAR file first!")
        
        self.poscar_info = {}

        with self.poscar_file.open() as f:
            self.poscar_info["title"] = f.readline().strip()
            if self.poscar_info["title"] == "":
                raise EOFError("Empty POSCAR file")
            
            try:
                self.poscar_info["scaling"] = float(f.readline().strip())
            except ValueError:
                raise ValueError("Unable to parse scaling factor")
            
            # check if there are 3 values
            # then check the data type of each
            try:
                self.poscar_info["lattice_x"] = list(map(float, f.readline().strip().split()))
                if len(self.poscar_info["lattice_x"]) < 3:
                    raise ValueError("Too few numbers supplied for lattice vector x")
                elif len(self.poscar_info["lattice_x"]) > 3:
                    raise ValueError("Too many numbers supplied for lattice vector x")
                
                for num in self.poscar_info["lattice_x"]:
                    float(num)
                
                self.poscar_info["lattice_y"] = list(map(float, f.readline().strip().split()))
                if len(self.poscar_info["lattice_y"]) < 3:
                    raise ValueError("Too few numbers supplied for lattice vector y")
                elif len(self.poscar_info["lattice_y"]) > 3:
                    raise ValueError("Too many numbers supplied for lattice vector y")
                
                for num in self.poscar_info["lattice_y"]:
                    float(num)
                
                self.poscar_info["lattice_z"] = list(map(float, f.readline().strip().split()))
                if len(self.poscar_info["lattice_z"]) < 3:
                    raise ValueError("Too few numbers supplied for lattice vector z")
                elif len(self.poscar_info["lattice_z"]) > 3:
                    raise ValueError("Too many numbers supplied for lattice vector z")
                
                for num in self.poscar_info["lattice_z"]:
                    float(num)
                
            except ValueError as err:
                raise ValueError("Unable to parse all lattice vectors: " + str(err))

            self.poscar_info["atoms"] = f.readline().strip().split()
            if len(self.poscar_info["atoms"]) < 1:
                raise ValueError("Unable to parse atom(s)")
            
            # can be two or more, corresponds with the number of unique atoms...
            try:
                self.poscar_info["atom_count"] = int(f.readline().strip())
            except ValueError:
                raise ValueError("Unable to parse atom count")
            
            self.poscar_info["coord_mode"] = f.readline().strip()

            # can TECHNICALLY use C, c, K, k, or anything else otherwise
            if self.poscar_info["coord_mode"] == "":
                raise ValueError("Unable to parse coordinate mode")
            elif (
                self.poscar_info["coord_mode"].lower() != "direct" and
                self.poscar_info["coord_mode"].lower() != "cartesian"
            ):
                raise ValueError("Invalid coordinate mode " + self.poscar_info["coord_mode"])
            
            for i in range(self.poscar_info["atom_count"]):
                atom_num = "atom_coord_" + str(i)
                self.poscar_info[atom_num] = list(map(float, f.readline().strip().split()))
                if len(self.poscar_info[atom_num]) < 3:
                    raise ValueError("Missing atom positions for atom " + str(i))
                elif len(self.poscar_info[atom_num]) > 3:
                    raise ValueError("Too many positions for atom " + str(i))
    
    def set_poscar(self, poscar_path: str) -> bool:
        """Set the POSCAR path
        Returns True if it is a valid POSCAR
        Return False if it is an invalid POSCAR"""
        self.file_validator = FileValidation(poscar_path)

        self.poscar_file = self.file_validator.valid_file()
        if self.poscar_file:
            try:
                self.__get_poscar_info()
            except ValueError as err:
                raise ValueError("Unable to parse POSCAR: " + str(err))
                
            print("Valid POSCAR!")
            return True
        else:
            print("Invalid POSCAR path given...")
            return False
        
    def set_potpawpbte(self, potpawpbte_path: str) -> bool:
        """Set the potpawpbte directory path
        Returns True if it is a valid POSCAR
        Return False if it is an invalid POSCAR"""
        self.file_validator = FileValidation(potpawpbte_path)

        self.potpawpbte_dir = self.file_validator.valid_directory()
        if self.potpawpbte_dir:
            return True
        else:
            print("Invalid directory path...")
            return False

    def valid_poscar_atoms(self) -> bool:
        """Validate that atoms in POSCAR are valid atoms and validate that the potpawpbte directory contains those atoms
        Returns True if all atoms are found in the potpawpbte directory
        Returns False if atoms could not be found in the potpawpbte directory"""
        success = False
        if self.poscar_info and self.potpawpbte_dir:
            found_atom = 0
            for atom in self.poscar_info["atoms"]:
                for child in self.potpawpbte_dir.iterdir(): 
                    split_path = str(child.absolute()).split("/")
                    atom_path = split_path[-1]
                    
                    if atom_path == atom:
                        found_atom += 1
                        if found_atom == len(self.poscar_info["atoms"]):
                            success = True
                            break
                else:
                    continue
                break
        return success


    def set_hdf5(self, hdf5_lib_path: str) -> bool:
        """Set the hdf5 library directory path
        Returns True if it is a valid directory
        Return False if it is an invalid directory"""
        self.file_validator = FileValidation(hdf5_lib_path)
        self.file_validator.change_path(hdf5_lib_path)

        self.hdf5_lib_dir = self.file_validator.valid_directory()

        if self.hdf5_lib_dir:
            print("Valid directory!\n")
            return True
        else:
            print("Invalid directory...\n")
            return False

    def set_incar(self, incar_path: str) -> bool:
        """Set the INCAR path
        Returns True if it is a valid file
        Return False if it is an invalid file"""
        incar_path = incar_path.strip()

        if not incar_path:
            print("Using default INCAR")

            incar_path = argv[0] + "/../default-input-files/INCAR"
        
        self.file_validator = FileValidation(incar_path)
        self.file_validator.change_path(incar_path)
        
        self.incar_file = self.file_validator.valid_file()
        if self.incar_file:
            print("Valid INCAR!\n")
            return True
        else:
            print("Invalid INCAR path given...\n")
            return False
    
    def set_vasp(self, vasp_bin_path: str) -> bool:
        """Set the VASP binary path
        Returns True if it is a valid file
        Return False if it is an invalid file"""
        self.file_validator.change_path(vasp_bin_path)

        self.vasp_bin_file = self.file_validator.valid_file()

        if self.vasp_bin_file:
            print("Valid VASP bin!\n")
            return True
        else:
            print("Invalid VASP bin...\n")
            return False
    
    def display_files(self):
        """Displays the files loaded into the class"""
        if not self.poscar_info:
            print("No POSCAR given, please initialize the class!")
            return
        print(
f"""
####################
    POSCAR INFO     
####################
Title: {self.poscar_info["title"]}

Lattice:
{self.poscar_info["lattice_x"]}
{self.poscar_info["lattice_y"]}
{self.poscar_info["lattice_z"]}
Lattice Scaling: {self.poscar_info["scaling"]}

Atoms: {self.poscar_info["atoms"]}
Atom Count: {self.poscar_info["atom_count"]}
Coordinate Mode: {self.poscar_info["coord_mode"]}
Atom Coordinates:
"""
        )
        for i in range(self.poscar_info["atom_count"]):
            atom_num = "atom_coord_" + str(i)
            print(self.poscar_info[atom_num])
    
    def __initialize(self):
        """REFACTORING IN PROGRESS"""

        init_dict = {}

        init_dict = self.__get_input_files()

        if init_dict != {}:

            # TODO: check if somehow the validation here fails. throw an exception
            print("Previous input files and paths detected, please review them below:")
            
            self.file_validator = FileValidation(init_dict["poscar"])
            self.poscar_file = self.file_validator.valid_file()

            self.file_validator.change_path(init_dict["incar"])
            self.incar_file = self.file_validator.valid_file()

            self.file_validator.change_path(init_dict["potpawpbte"])
            self.potpawpbte_dir = self.file_validator.valid_directory()

            self.file_validator.change_path(init_dict["hdf5-lib"])
            self.hdf5_lib_dir = self.file_validator.valid_directory()

            self.file_validator.change_path(init_dict["vasp_std"])
            self.vasp_bin_file = self.file_validator.valid_file()

            self.__get_poscar_info()

            print("POSCAR: ", init_dict["poscar"])
            print("INCAR: ", init_dict["incar"])
            print("potpawpbte: ", init_dict["potpawpbte"])
            print("hdf5 library: ", init_dict["hdf5-lib"])
            print("vasp_std bin: ", init_dict["vasp_std"])
            return
        
        init_dict = self.__get_input_files()

        print("Please review the input files below:")
        print("POSCAR: ", init_dict["poscar"])
        print("INCAR: ", init_dict["incar"])
        print("potpawpbte: ", init_dict["potpawpbte"])
        print("hdf5 library: ", init_dict["hdf5-lib"])
        print("vasp_std bin: ",init_dict["vasp_std"])
        
    def __write_out_files(self):
        """REFACTORING IN PROGRESS"""
        init_json = Path(argv[0] + "/../default-input-files/init.json").resolve()
        with init_json.open('w') as f:
            json.dump(
                {
                    "poscar": str(self.poscar_file),
                    "potpawpbte": str(self.potpawpbte_dir),
                    "hdf5-lib": str(self.hdf5_lib_dir),
                    "incar": str(self.incar_file),
                    "vasp_std": str(self.vasp_bin_file)
                }, 
                f, sort_keys = True, indent = 4
            )
    
class VaspRelax(VaspInit):
    """Handles relaxation calculations"""

    def perform_calculation(self, calc_dir: str):
        # create relaxation directory

        # copy POSCAR, INCAR and create POTCAR in that directory
        result = subprocess.run(["cp", str(self.poscar_file), calc_dir + "/POSCAR"], capture_output = True, text = True)
        result = subprocess.run(["cp", str(self.incar_file), calc_dir + "/INCAR"], capture_output = True, text = True)

        new_potcar = Path(str(calc_dir) + "/POTCAR")

        new_potcar.touch()

        with new_potcar.open('w') as f:
            cat_command = ["cat"]
            for atom in self.poscar_info["atoms"]:
                cat_command.append(str(self.potpawpbte_dir) + "/" + atom + "/POTCAR")
                
            subprocess.run(cat_command, stdout = f)
        
        # create script file
        relaxation_script = Path(calc_dir + "/job_vasp66_coeus")

        relaxation_script.touch(0o755)

        # can definitely prompt for info in this
        relaxation_script.write_text(
f"""#!/bin/bash
#SBATCH --partition short      ### Partition (short/60min, medium/4days/dual10, long/20days/dual10, himem/20days/dual12, phi/20days/64cores)
#SBATCH --job-name relaxation  ### Job Name
#SBATCH --nodes 2              ### Number of Nodes (14 cores/cpu, 2 cpu/node. Remember NPAR=14)
#SBATCH --ntasks-per-node 20   ### Number of OpenMP threads (28 maximum, 14 for underpopulation to get more memory
#SBATCH --time 01:00:00
#SBATCH --output output        ### Note that cpus-per-task can be any even number but you must then adjust NPAR to be a factor of MPIPROCS)
#SBATCH --error vasp-test-%j.err
#SBATCH --output vasp-test-%j.out

module load oneapi/2024.2.0
export LD_LIBRARY_PATH={str(self.hdf5_lib_dir)}/build/lib:$LD_LIBRARY_PATH

ulimit -s unlimited
vasp={str(self.vasp_bin_file)}

mpirun -n 40 $vasp > stdout-vasp  ### Standard VASP compilation = _std, NonCollinear = _ncl
""")

        script = Path(argv[0] + "/../script-files/start_relaxation").resolve()
        subprocess.run(["bash", str(script), calc_dir])

class VaspPhelel(VaspInit):
    """Handles electron-phonon interaction calculations"""

    def perform_calculation(self, relaxed_dir: str, calc_dir: str):
        # move POSCAR and POTCAR to the new directory
        subprocess.run(["cp", relaxed_dir + "/CONTCAR", calc_dir + "/POSCAR"])
        subprocess.run(["cp", relaxed_dir + "/POTCAR", calc_dir + "/POTCAR"])

        template_toml = Path(calc_dir + "/template.toml")

        template_toml.write_text(f"""# Template for velph
 
[init.options]
kspacing = 0.2
kspacing_dense = 0.14
tolerance = 1e-4

# supercell_matrix is dynamic per POSCAR
# max number of atoms is 60 but max_num_atoms may not
# product the desired supercell dimension 
# (would need further testing)

[vasp.incar]
encut = 520
kspacing = 0.2
ncore = 2

[vasp.phelel.incar]
elph_prepare = true
lwap = false
kspacing = 0.2
gga = "PS"

[vasp.transport.incar]
elph_selfen_temps = 300
gga = "PS"

[scheduler]
p_code = "short"
job_name = "vasp-elph"
mpirun_command = "mpirun -n 40"
vasp_binary = "$vasp > stdout-vasp"
modules = \"\"\"\\
#module purge
#module load openmpi-2.0/gcc-7.2.0
#module load Chemistry/VASP/5.4.4/openmpi-2.0/gcc-7.2.0
\"\"\"
prepend_text = \"\"\"\\
module load oneapi/2024.2.0
export LD_LIBRARY_PATH={str(self.hdf5_lib_dir)}/build/lib:$LD_LIBRARY_PATH

ulimit -s unlimited
vasp={str(self.vasp_bin_file)}
\"\"\"
append_text = \"\"
nodes = 2
ntasks = 20
walltime = "4:00:00"
scheduler_template = \"\"\"\\
#!/bin/bash
#SBATCH --partition={{p_code}}
#SBATCH --job-name={{job_name}}
#SBATCH --nodes={{nodes}}
#SBATCH --ntasks-per-node={{ntasks}}
#SBATCH --time={{walltime}}
#SBATCH --error vasp-test-%j.err
#SBATCH --output vasp-test-%j.out

{{modules}}

{{prepend_text}}
{{mpirun_command}} {{vasp_binary}} > stdout-vasp
{{append_text}}
\"\"\"

        """)
        script = Path(argv[0] + "/../script-files/start_phelel").resolve()
        subprocess.run(["bash", str(script), calc_dir])


class VaspTransport(VaspInit):
    """Handles transport property calculations"""

    def perform_calculation(self, phelel_dir: str):
        # given a directory with completed dispacement directories, do the transport property stuff
        script = Path(argv[0] + "/../script-files/start_transport").resolve()
        subprocess.run(["bash", str(script), phelel_dir])

class VaspInterface():
    """Interface class that interacts with the hierarchy and handles multiple relaxations"""
    def __init__(self):
        self.vasp_relax = VaspRelax() 
        self.vasp_phelel = VaspPhelel()
        self.vasp_transport = VaspTransport()
    
    def __get_first_word_last_line(self, calc_dir: str):
        """Extract the first word from the last line of the OUTCAR"""
        
        outcar_path = Path(calc_dir + "/OUTCAR")
        if not outcar_path.exists():
            raise SystemError("OUTCAR in the given directory does not exist")
        
        return subprocess.run(["tail", "-n", "1", str(outcar_path)], capture_output = True, text = True)
    
    
    def __determine_continuation(self, log_dictionary: dict, step_name: str) -> bool:
        """Determines if the program should continue onto the next step in calculations
        Returns True if there are successful calculations
        Returns False if there are no more successful calculations"""
        log_poscars = list(log_dictionary.keys())
        for file in log_poscars:
            if log_dictionary[file] == "F":
                del log_dictionary[file]
        
        # check if the dictionary is empty, return!
        if log_dictionary == {}:
            print(f"ERROR: No jobs were successful during the {step_name} step! Exiting...")
            return False
        return True
    
    def __log_in_progress_directory(self, log_dictionary: dict, log: Path):
        """Write out and \"L\" for loading and processing the calculation directory."""
        for poscar_path in log_dictionary.keys():
            
            # print each of the children to the log file
            with log.open('a') as f:
                f.write("L " + str(poscar_path) + "\n")
    
    def __create_new_calculation_path(self, calc_path: str, poscar_path: str) -> Path:
        """Creates a new Path given the base calculation path and the path with the POSCAR name"""
        # return Path("calc_path/POSCAR_NAME")
        return Path(str(calc_path) + "/" + str(poscar_path).split('/')[-1])
    
    def __log_calculation_status(self, calc_dir: Path) -> str:
        """Output the status of the calculation directory
        Returns \"P\" if the calculation passed
        Returns \"F\" if the calculation failed"""
        outcar_path = Path(calc_dir) / "OUTCAR"

        try:
            result = self.__get_first_word_last_line(str(calc_dir))
            last_word = result.stdout.strip().split()
            try:
                last_word = last_word[0]
            except:
                # the last line was empty, meaning it didn't finish in time
                return "F"
            
            # calculation finished woohoo
            if last_word == "Voluntary":
                return "P"
            else:
                # OUTCAR didn't finish in time
                return "F"
        except:
            # OUTCAR was never made
            return "F"
    
    # TODO: break down this function and rename it to be more accurate to what it does
    def initialize(self):
        """The main application function"""
        valid_poscar_directory: bool = False
        # TODO: make this a function
        while not valid_poscar_directory:
            poscar_dir = input("Enter the full directory path with the POSCARs: ")

            file_validator = FileValidation(poscar_dir)
            valid_poscar_directory = file_validator.valid_directory()
        
        poscar_dir = Path(poscar_dir)

        valid_potpawpbte: bool = False
        while not valid_potpawpbte:
            potpawpbte_path = input("Enter the directory that holds the potpawpbte: ")
            valid_potpawpbte = self.vasp_relax.set_potpawpbte(potpawpbte_path)

        valid_hdf5: bool = False
        while not valid_hdf5:
            hdf5_lib_path = input("Enter the hdf5-1.14.6 directory path: ")
            valid_hdf5 = self.vasp_relax.set_hdf5(hdf5_lib_path)

        valid_incar: bool = False
        while not valid_incar:
            # now create relaxation INCAR
            incar_path = input("Enter the INCAR file path to use (leave blank to use default): ")
            valid_incar = self.vasp_relax.set_incar(incar_path)

        valid_vasp: bool = False
        while not valid_vasp:
            vasp_path = input("Enter the vasp bin path to use: ")
            valid_vasp = self.vasp_relax.set_vasp(vasp_path)
        
        # ask for max time in minutes for user experience
        relax_time = int(input("Enter the max time in minutes for the relaxation step: "))
        phelel_time = int(input("Enter the max time in minutes for the phelel step: "))
        transport_time = int(input("Enter the max time in minutes for the transport step: "))

        # convert all of them into seconds for consistency with time module
        relax_time *= 60
        phelel_time *= 60
        transport_time *= 60

        print("\nThe following folders will be created:")
        print(Path(str(poscar_dir) + "/../poscar_relaxed").resolve())
        print(Path(str(poscar_dir) + "/../poscar_phelel").resolve())
        print("Ensure these directories do not exist before continuing.")
        input("Press ENTER to continue: ")

        print("Creating /poscar_relaxed folder as sibling to supplied folder...")
        poscar_relax_dir = Path(str(poscar_dir) + "/../poscar_relaxed").resolve()
        poscar_relax_dir.mkdir()
        poscar_relax_log = Path(str(poscar_dir) + "/../poscar_relaxed/relaxed.log").resolve()
        original_poscar_dir = sorted(poscar_dir.iterdir())
        
        # set up the dictionary holding the poscar and status of the calculation
        log_dictionary = {}
        
        for poscar_path in original_poscar_dir:
            log_dictionary[poscar_path] = "L"

        self.__log_in_progress_directory(log_dictionary, poscar_relax_log)

        # NOTE: this part can probably be split into either processes or threads
        for poscar_path in original_poscar_dir:
            # feed in a directory for it to copy everything into
            calc_dir = self.__create_new_calculation_path(poscar_relax_dir, poscar_path)
            
            calc_dir.mkdir()
            try:
                valid_poscar = self.vasp_relax.set_poscar(str(poscar_path))
            except Exception as err:
                print(f"WARNING: {poscar_path}: {err}")
            else:
                if self.vasp_relax.valid_poscar_atoms():
                    self.vasp_relax.perform_calculation(str(calc_dir))
        check_script = Path(argv[0] + "/../script-files/check_queue").resolve()

        start_time = time.time()
        result = subprocess.run(["bash", str(check_script)], capture_output = True, text = True)

        try:
            job_count = int(result.stdout) - 1
        except:
            print(f"ERROR: No valid POSCARS were submitted for a job. Exiting...")
            return
        else:
            current_time = time.time()
            time_delta = current_time - start_time
            while job_count >= 1 and time_delta < relax_time:
                time.sleep(300)
                current_time = time.time()
                result = subprocess.run(["bash", str(check_script)], capture_output = True, text = True)
                job_count = int(result.stdout) - 1
            
        # TODO: once this is working see if you can just open the file once
        with poscar_relax_log.open('w') as f:
            f.write("")
        
        # -------- after relaxation step -------
        # create the dictionary
        log_dictionary = {}

        # TODO: the program needs to check the OUTCAR for the pullay force and external pressure.
        # if its too big then the calculation needs to be reran.
        # may need a grep to see the external pressure and pullay to see if it needs more calculation
        # AFTER I GET THE OTHER STEPS WORKING

        # check the state of all the programs. They need to have an OUTCAR that has the time completed near the end.
        for calc_dir in sorted(poscar_relax_dir.iterdir()):
            if calc_dir.is_file():
                continue
            
            relax_status = self.__log_calculation_status(calc_dir)
            with poscar_relax_log.open('a') as f:
                f.write(relax_status + " " + str(calc_dir) + "\n")
            log_dictionary[calc_dir] = relax_status

        if self.__determine_continuation(log_dictionary, "relaxation") == False:
            return
        
        # set up the directory. the root directory for which all calculations occur will be given by the user
        poscar_phelel_dir = Path(str(poscar_dir) + "/../poscar_phelel").resolve()
        poscar_phelel_dir.mkdir()
        poscar_phelel_log = Path(str(poscar_dir) + "/../poscar_phelel/phelel.log").resolve()

        self.__log_in_progress_directory(log_dictionary, poscar_phelel_log)
        
        # copy over the relax parameters to the phelel object
        self.vasp_phelel.copy_init(self.vasp_relax)

        for relaxed_directory in log_dictionary.keys():
            # use the name of the original poscar as the directory name
            calc_dir = self.__create_new_calculation_path(poscar_phelel_dir, relaxed_directory)
            
            # make new directories based off the ones that were successful.
            calc_dir.mkdir()
            
            self.vasp_phelel.perform_calculation(str(relaxed_directory), str(calc_dir))

        start_time = time.time()
        result = subprocess.run(["bash", str(check_script)], capture_output = True, text = True)
        job_count = int(result.stdout) - 1
        current_time = time.time()
        time_delta = current_time - start_time
        while job_count >= 1 and time_delta < phelel_time:
            time.sleep(300)
            current_time = time.time()
            result = subprocess.run(["bash", str(check_script)], capture_output = True, text = True)
            job_count = int(result.stdout) - 1
        
        # reset the file
        with poscar_phelel_log.open('w') as f:
            f.write("")
        
        # Refresh the dictionary
        log_dictionary = {}

        # iterate through the directories in poscar_phelel
        # This is probably n x 3 so 3n or n
        # for each folder in poscar_phelel
        for sub_dir in sorted(poscar_phelel_dir.iterdir()):
            if sub_dir.is_file():
                continue
            
            phelel_status = "P"
            phelel_calc_dir = Path(str(sub_dir) + "/phelel")
            for disp_dir in sorted(phelel_calc_dir.iterdir()):
                if disp_dir.is_file():
                    continue
                
                # check if there is an OUTCAR for each disp-### directory
                phelel_status = self.__log_calculation_status(disp_dir)
                
                if phelel_status == "F":
                    break

            with poscar_phelel_log.open('a') as f:
                f.write(phelel_status + " " + str(sub_dir) + "\n")
            log_dictionary[sub_dir] = phelel_status
                
        if self.__determine_continuation(log_dictionary, "phelel") == False:
            return
        
        # ---- transport calculation step ---
        # now its time for the transport calculation!
        vasp_transport.copy_init(self.vasp_phelel)
        
        # it will reuse a directory, but the transport property log will be in there
        # this is to reuse directories and not increase space unnecessarily

        poscar_transport_log = Path(str(phelel_dir) + "/../poscar_phelel/transport.log").resolve()
        
        self.__log_in_progress_directory(log_dictionary, poscar_transport_log)

        for phelel_directory in log_dictionary.keys():
            self.vasp_transport.perform_calculation(str(phelel_directory))


        start_time = time.time()
        result = subprocess.run(["bash", str(check_script)], capture_output = True, text = True)
        job_count = int(result.stdout) - 1
        current_time = time.time()
        time_delta = current_time - start_time
        while job_count >= 1 and time_delta < transport_time:
            time.sleep(300)
            current_time = time.time()
            result = subprocess.run(["bash", str(check_script)], capture_output = True, text = True)
            job_count = int(result.stdout) - 1

        # reset the file
        with poscar_transport_log.open('w') as f:
            f.write("")
        
        # the program has to go by the directories in the transport log and can't just iterate
        # over the folders as usual because some of them failed before

        for phelel_dir in log_dictionary.keys():
            transport_dir = Path(str(phelel_dir) + "/transport")
            
            transport_status = self.__log_calculation_status(transport_dir)
            with poscar_transport_log.open('a') as f:
                f.write(transport_status + " " + str(phelel_dir) + "\n")
            log_dictionary[phelel_dir] = transport_status

        print("All calculations complete; please check the log files.")
