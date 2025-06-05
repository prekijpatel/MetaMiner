"""
MetaMiner - Streamlined GUI Tool for Retrieving, Normalizing and Exploring Metadata
Copyright (C) [2025] [Patel Jaykumar Kiritkumar, 
Molecular Imaging and Molecular Diagnostics Lab, 
Indian Institute of Technology Delhi, India]

Licensed under the GNU Affero General Public License v3.0 (AGPL-3.0)
with additional restrictions: **Non-Commercial Use Only**.

For full terms, see the LICENSE and LICENSE-EXCEPTIONS.md files.
For commercial licensing inquiries, contact: prekijpatel2.0@gmail.com
"""

import os
import multiprocessing
import sys
import pickle
import threading
from typing import Any, Callable
from tkinter import messagebox, filedialog
import customtkinter
from utils import *
import logging

sys.stdout = sys.stderr = open(os.devnull, 'w')

logging.basicConfig(filename="log_file.log", level=logging.INFO, format='%(asctime)s:%(levelname)s - %(message)s')

class gui_button(customtkinter.CTkButton):
    def __init__(self, master, text: str, command: Callable[[], Any], height: int, width: int, corner_radius: int, text_color: str, fg_color: str, border_color: str, border_width: int):
        
        font = customtkinter.CTkFont(size=12, weight="bold")
        
        super().__init__(master, text=text, command=command, height=height, width=width, corner_radius=corner_radius, text_color=text_color, fg_color=fg_color, border_color=border_color, border_width=border_width, font=font)

class gui_frame_labels(customtkinter.CTkLabel):
    def __init__(self, master, text: str, fg_color: str, text_color: str, compound: str = 'left', padx: int = 5, pady: int = 5, bold: bool = False):

        font = customtkinter.CTkFont(size=14, weight="bold") if bold else customtkinter.CTkFont(size=14)
        
        super().__init__(master, text=text, fg_color=fg_color, text_color=text_color, font=font, compound=compound, padx=padx, pady=pady)

class gui_dropdown(customtkinter.CTkOptionMenu):
    def __init__(self, master, values, command: Callable[[], Any]):
        super().__init__(master, values=values, command=command, width=150, height=30, corner_radius=6, anchor='center', text_color='#FFFFFF', fg_color='#000000', button_color='#000000', button_hover_color='#5C2751', dropdown_fg_color='#000000', dropdown_hover_color="#5C2751", dropdown_text_color='#FFFFFF')

class gui_title(customtkinter.CTkLabel):
    def __init__(self, master, text: str, fg_color: str, text_color: str, compound: str = 'left', padx: int = 5, pady: int = 5):

        font = customtkinter.CTkFont(size=16, weight="bold")
        
        super().__init__(master, text=text, fg_color=fg_color, text_color=text_color, font=font, compound=compound, padx=padx, pady=pady)

class gui_labels(customtkinter.CTkLabel):
    def __init__(self, master, text: str, fg_color: str, text_color: str, compound: str = 'left', padx: int = 5, pady: int = 5, bold: bool = False, corner_radius: int = 0, image: str = None, width: int = 0, height: int = 0, wraplength: int = 0):

        font = customtkinter.CTkFont(size=12, weight="bold") if bold else customtkinter.CTkFont(size=12)
        
        super().__init__(master, text=text, fg_color=fg_color, text_color=text_color, font=font, compound=compound, padx=padx, pady=pady, corner_radius=corner_radius, image=image, width=width, height=height, wraplength=wraplength)


class gui_textbox(customtkinter.CTkTextbox):
    def __init__(self, master, width: int, height: int, padx: int, pady: int, corner_radius: int, text_color: str, fg_color: str, border_color: str, border_width: int, wrap: str = "none"):
        super().__init__(master, width=width, height=height, padx=padx, pady=pady, corner_radius=corner_radius, text_color=text_color, fg_color=fg_color, border_color=border_color, border_width=border_width, wrap=wrap)


class logging_window(customtkinter.CTkToplevel):
    def __init__(self, master, title: str, width: int, height: int, bg_color: str, fg_color: str):
        super().__init__(master)
        self.title(title)
        self.geometry(f"{width}x{height}")
        self.configure(bg=bg_color)
        self.fg_color = fg_color
        self.iconbitmap(os.path.abspath('./img/meta.ico'))
        

class LogHandler(logging.Handler):
    def __init__(self, text_widget):
        super().__init__()
        self.text_widget = text_widget

    def emit(self, record):
        msg = self.format(record)
        # Append the message to the text widget
        self.text_widget.insert(customtkinter.END, f"{msg}\n")
        self.text_widget.see(customtkinter.END)


class metaminer(customtkinter.CTk):
    def __init__(self):
        super().__init__(fg_color=('black'))
        # for light theme
        # super().__init__(fg_color=('white'))
        
        self.geometry('600x600')
        self.maxsize(600, 600)
        self.minsize(600, 600)
        self.title('MetaMiner')

        self.iconbitmap(os.path.abspath('./img/meta.ico'))

        logging.debug("Configuring the window...")

        self.log = threading.Thread(target=self.open_log_window)
        self.log.start()

        self.save_path = threading.Thread(target=self.ask_for_file_path)
        self.save_path.start()

        
        # clean up method incorporated
        # self.protocol("WM_DELETE_WINDOW", self.cleanup)

        # configuring entire window into two spaces
        self.rowconfigure((0,1,2,3), weight=1)
        self.rowconfigure(4, weight=7)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        
        # <<< Top shelf of the window starts here. >>>
        # In top shelf of the window, putting the load JSON button.
        
        self.welcome_text = gui_title(
            self, 
            text="Welcome to Metaminer v0.1",
            fg_color='black',
            text_color='#8EA8FF'
        )
        
        self.welcome_text.grid(
            row=0, column=0, columnspan=2, padx = 7, pady = 7
        )
        
        self.load_json_text = gui_frame_labels(
            self, 
            text="Choose metadata from computer:", 
            fg_color='black', 
            text_color='white', 
            bold=True
        )    
        self.load_json_text.grid(
            row=1, column=0, padx = 7, pady = 7, sticky="e"
        )

        # text box for showing/typing the path of the selected file
        self.textbox = gui_textbox(
            self, 
            width=200, 
            height=5, 
            corner_radius=6, 
            text_color='white', 
            fg_color='#141414', 
            border_color='grey', 
            border_width=1, 
            padx=5, 
            pady=5
        )
        self.textbox.grid(
            row=1, column=1, padx = 1, pady = 1, sticky="w"
        )

        # button for choosing the file
        self.choose_file_button = gui_button(
            self, 
            text='Choose File', 
            command=self.choose_file, 
            height=15, 
            width=100, 
            corner_radius=5, 
            text_color='white', 
            fg_color='#141414', 
            border_color='skyblue', 
            border_width=1
        )
        self.choose_file_button.grid(
            row=2, column=1, columnspan=2, padx = (35, 10), pady = 1, sticky="n", ipadx=5
        )

        # button for starting the analysis for the selected file
        self.analyze_button_top = gui_button(
            self, 
            text='Analyze... { }', 
            command=self.analyze_button_top, 
            height=28, 
            width=140, 
            corner_radius=6, 
            text_color='black', 
            fg_color='#8EA8FF', 
            border_color='#141414', 
            border_width=1
        )
        self.analyze_button_top.grid(row=3, column=0, columnspan=2, padx = 1, pady = 1)
        # <<< Top shelf of the window ends here. >>>

        # <<< Bottom shelf of the window starts here. >>>
        # Now, in bottom shelf of the window, putting various dropdowns for downloading json based on user input.
        
        ## creating frame for better look and separation of upper and lower shelf
        self.download_json_frame = customtkinter.CTkFrame(self, fg_color="#141414")
        
        self.download_json_frame.grid(
            row=4, column=0, columnspan=2, padx = 7, pady = 7, sticky='nsew'
        )

        ## creating dropdown for selecting organism
        
        ### spacing things for better look - very neive, I know!
        self.download_json_frame.rowconfigure(
            (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30), weight=1
        )

        self.download_json_label = gui_frame_labels(
            self.download_json_frame, text="Load metadata from NCBI server:", text_color='white', fg_color='#141414', bold=True
        )
        self.download_json_label.grid(
            row=0, column=0, columnspan=3, padx = 10, pady = 6, sticky="w"
        )

        self.download_json_frame.columnconfigure(0, weight=1)
        self.download_json_frame.columnconfigure(1, weight=2)
        self.download_json_frame.columnconfigure(2, weight=2)
        self.download_json_frame.columnconfigure(3, weight=2)
        self.download_json_frame.columnconfigure(4, weight=2)
        self.download_json_frame.columnconfigure(5, weight=1)

        # ### Initilizing `information` clipart to go with the dropdown labels

        # info_clipart = customtkinter.CTkImage(
        #     Image.open('./img/info_clipart.png'), size=(20, 20)
        # )

        # #### creating tooltips for Organism, Genome/Gene and Accession/Taxon dropdown menus
        # self.Organism_tooltip = gui_labels(
        #     self.download_json_frame, 
        #     text="",
        #     fg_color="#FFFFC5", 
        #     text_color='black',
        #     corner_radius=5,
        #     width=200,
        #     height=50,
        #     wraplength=40
        # )
        
        # self.genome_gene_tooltip = gui_labels(
        #     self.download_json_frame, 
        #     text="",
        #     fg_color="#FFFFC5", 
        #     text_color='black',
        #     corner_radius=5,
        # )

        # self.accid_taxid_tooltip = gui_labels(
        #     self.download_json_frame, 
        #     text="",
        #     fg_color="#FFFFC5", 
        #     text_color='black',
        #     corner_radius=5,
        # )

        ### creating label and dropdown for selecting organism
        self.Organism_label = gui_labels(
            self.download_json_frame,
            text="Organism",
            text_color='white',
            fg_color='#141414',
            bold=True,
            # image=info_clipart,
        )
        self.Organism_label.grid(row=1, column=1, padx = 10, pady = 6)
        # self.Organism_label.bind("<Enter>", self.show_organism_tooltip)
        # self.Organism_label.bind("<Leave>", self.hide_organism_tooltip)
        
        self.Organism_menu1 = gui_dropdown(
            self.download_json_frame, 
            values=['--select--','prokaryote', 'eukaryote', 'virus'], 
            command=self.get_organism
        )
        self.Organism_menu1.set('--select--')
        self.Organism_menu1.grid(
            row=2, column=1, padx = 10, pady = 0
        )

        ### creating label and dropdown for selecting genome/gene

        self.genome_gene_label = gui_labels(
            self.download_json_frame,
            text="Genome/Gene",
            text_color='white',
            fg_color='#141414',
            bold=True,
        )
        self.genome_gene_label.grid(
            row=1, column=2, padx = 10, pady = 6
        )

        self.genome_gene_menu2 = gui_dropdown(
            self.download_json_frame,
            values=[],
            command=self.get_genome_gene
        )
        self.genome_gene_menu2.set('--select--')
        self.genome_gene_menu2.grid(
            row=2, column=2, padx = 10, pady = 0
        )

        ### creating label and dropdown for selecting accession/taxon

        self.accid_taxid_label = gui_labels(
            self.download_json_frame,
            text="Accession/Taxon",
            text_color='white',
            fg_color='#141414',
            bold=True,    
        )
        self.accid_taxid_label.grid(
            row=1, column=3, padx = 10, pady = 6
        )

        self.accid_taxid_menu3 = gui_dropdown(
            self.download_json_frame,
            values=[],
            command=self.get_accid_taxid
        )
        self.accid_taxid_menu3.set('--select--')
        self.accid_taxid_menu3.grid(
            row=2, column=3, padx = 10, pady = 0
        )

        ### now putting text box to enter the accession number
        self.accession_label = gui_labels(
            self.download_json_frame,
            text="Enter details:",
            text_color='white',
            fg_color='#141414',
            bold=True
        )
        self.accession_label.grid(
            row=1, column=4, padx = 10, pady = 6
        )

        self.inputtextbox = gui_textbox(
            self.download_json_frame, width=140, height=4, corner_radius=6, text_color='white', fg_color='black', border_color='black', border_width=1, padx=5, pady=5
        )
            
        self.inputtextbox.grid(
            row=2, column=4, padx = 10, pady = 0
        )

        ## now putting the download button

        self.download_json_button = gui_button(
            self.download_json_frame, 
            text='Download JSON { }', 
            command=self.download_json_button, 
            height=28, 
            width=140, 
            corner_radius=6, 
            text_color='black', 
            fg_color='#8EA8FF', 
            border_color='#141414', 
            border_width=1
        )
        self.download_json_button.grid(
            row=8, column=0, columnspan=5, padx = 1, pady = 1
        )
        
        self.analyze_button_bottom = gui_button(
            self.download_json_frame, 
            text='Analyze... { }', 
            command=self.analyze_button_bottom, 
            height=28, 
            width=140, 
            corner_radius=6, 
            text_color='black', 
            fg_color='#8EA8FF', 
            border_color='#141414', 
            border_width=1
        )
        self.analyze_button_bottom.grid(
            row=9, column=0, columnspan=5, padx = 1, pady = 1
        )

    
        # <<< Bottom shelf of the window ends here. >>>

        # <<< Setting up the debug options. >>>
        # --- Debug checkbox in bottom-right corner ---
        self.debug_var = customtkinter.StringVar(value="off")
        self.debug_checkbox = customtkinter.CTkCheckBox(
            self.download_json_frame,
            text="Debug Mode",
            variable=self.debug_var,
            onvalue="on",
            offvalue="off",
            command=self.toggle_debug_mode,
            checkbox_height=15,
            checkbox_width=15,
            border_width=1.5,
            corner_radius=1,
            hover_color="#5C2751",
            text_color="white",
        )
        self.debug_checkbox.place(relx=1.0, rely=1.0, anchor="se", x=-10, y=-10) 





    # <<< General functions to be used in the window. >>>

    def perform_analysis(self):
        self.count = record_numbers(self.data)
        # self.biosample_att = biosample_attributes_from_json(self.data, self.count)
        # self.assembly_att = assembly_attributes_from_json(self.data, self.count)
        # self.annotation_att = annotation_attributes_from_json(self.data, self.count)
        self.all_att = get_allmetadata(self.data, self.count, self.saving_file_path)
    
    def load_and_analyze(self):
        try:
            self.data = load_json(self.file_path)
            
            self.result = self.perform_analysis()

            if self.all_att is not None and not self.all_att.empty:
                self.good_for_view = True
                logging.debug(f"good_for_view set to {self.good_for_view}. Loding explore_data button.")

                if self.good_for_view:
                    self.explore_button = gui_button(
                        self.log_window, text="Explore Data", command=self.explore_data,
                        height=40, width=150, corner_radius=8, text_color="white",
                        fg_color="#5C2751", border_color="#5C2751", border_width=2
                    )
                    self.explore_button.pack(pady=10)
            else:
                self.good_for_view = False

        except Exception as e:
            logging.error(f"Error during analysis: {e}")
            self.after(0, lambda: messagebox.showerror("Error", "An error occurred!"))
    
    def explore_data_commands(self):
        # self.logger.addHandler(self.log_handler)

        logging.info("Starting normalization and categorization of data for visualization and exploration.")

        # Ensure `self.all_att` is a DataFrame
        if isinstance(self.all_att, pd.DataFrame):
            try:
                # Normalize and categorize the data
                self.df = all_normalization_operations(self.all_att, self.saving_file_path)
                
                # Ensure `self.df` is also a DataFrame
                if isinstance(self.df, pd.DataFrame):
                    logging.info("Normalization and categorization complete!")
                    try:
                        # # logging.info("Starting the Dash app!")

                        # self.start_dash()

                        # # self.run_dash = threading.Thread(target=meta_mined, args=(self.df,))
                        # # self.run_dash = multiprocessing.Process(target=meta_mined, args=(self.df,))
                        # # self.run_dash.start()    
                        # # meta_mined(self.df)                    

                        self.run_dash = threading.Thread(target=self.start_dash)
                        self.run_dash.start()


                    except Exception as e:
                        logging.error(f"Error during categorization: {e}", exc_info=True)
                        messagebox.showerror("Error", f"Error during categorization: {e}")
                else:
                    logging.error("Dataframe not found!")
                    messagebox.showerror("Error", "Dataframe not found!")
            except Exception as e:
                logging.error(f"Error during normalization: {e}", exc_info=True)
                messagebox.showerror("Error", f"Error during normalization: {e}")
        else:
            logging.error("Data not found!")
            messagebox.showerror("Error", "Data not found!")

    def start_dash(self):
        logging.info("Starting the dash app!")
        # try:
        #     # You can use multiprocessing or threading based on your requirement
        #     if threading.current_thread() == threading.main_thread():
        #         main_thread = "True"
        #     else:
        #         main_thread = "False"
            
        #     logging.debug(f"Starting the dashboard in main thread: {main_thread}")
        #     logging.debug(f"main thread:{threading.main_thread()}")
        #     # logging.debug(f"current process:{multiprocessing.current_process()}")
        #     logging.debug(f"current thread:{threading.current_thread()}")

        #     with open("tmp/df.pkl", "wb") as f:
        #         pickle.dump(self.df, f)
            
        #     with open("tmp/save_path.pkl", "wb") as f:
        #         pickle.dump(self.saving_file_path, f)

        #     self.run_dash_app()
        #     # self.run_dash = multiprocessing.Process(target=self.run_dash_app)
        #     # self.run_dash.start()
        
        # except Exception as e:
        #     logging.error(f"Error launching dashboard: {e}", exc_info=True)
        #     messagebox.showerror("Error", f"Error launching dashboard: {e}")
        
        try:
            logging.debug("We are here with dataframe and loading it into meta_mined function.")
            meta_mined_app(self.df, self.saving_file_path)
        except Exception as e:
            logging.debug(f"An error ocuured while running the dash app: {e}")

        #     meta_mined(self.df, self.saving_file_path)
    # def run_dash_app(self):
        
        # subprocess.run([sys.executable, "dashboard.py"], env={**os.environ, "PYINSTALLER_RESET_ENVIRONMENT": "1"})
        # subprocess.Popen([sys.executable], env={**os.environ, "PYINSTALLER_RESET_ENVIRONMENT": "1"})
        # sys.exit(0)
    
    def explore_data(self):
        logging.debug("explore_data button pressed!")
        threading.Thread(target=self.explore_data_commands).start()

    def update_label(self, text):
        self.process_label.configure(text=text)
    
    def open_log_window(self):
        # initializing the window for logging
        self.log_window = logging_window(
            self,
            title="metaminer is mining!",
            width=800,
            height=500,
            bg_color="#1A1A1A", 
            fg_color="#000000",
        )

        # putting a textbox in it
        self.log_textbox = gui_textbox(
            self.log_window, width=780, height=400, padx=10, pady=10,
            corner_radius=8, text_color="white", fg_color="#333333",
            border_color="#5C2751", border_width=2, wrap="word"
        )
        self.log_textbox.pack(expand=True, fill="both", padx=10, pady=10)

        # setting up the logger
        self.log_handler = LogHandler(self.log_textbox)
        self.log_handler.setFormatter(logging.Formatter('%(asctime)s:%(levelname)s - %(message)s'))
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(self.log_handler)
        
    def ask_for_file_path(self):
        try:
            self.saving_file_path = filedialog.askdirectory()
            logging.info(f"Saving file path selected: {self.saving_file_path}")
        except:
            self.saving_file_path = None
    
    def cleanup(self):
        # logging.info("Stopping MetaMiner...")
        try:
            if self.run_dash and self.run_dash.is_alive():
                # logging.info("Stopping the Dash server...")
                try:
                    if meta_mined_app:
                        meta_mined_app.shutdown()  
                    # logging.info("Dash server stopped.")
                except Exception as e:
                    # logging.error(f"Error while shutting down Dash server: {e}")
                    pass
                try:
                    if self.run_dash:
                        self.run_dash.join(timeout=5)
                    else:
                        pass
                except Exception as e:
                    pass
                    # logging.error(f"Error while joining Dash server thread: {e}")
                #     pass        
                # if self.run_dash.is_alive():
                #     # logging.warning("Dash server did not stop within the timeout.")
                #     pass
                # else:
                #     # logging.info("Dash server stopped.")
                #     pass
        except:
            pass
        
        try:
            self.quit()
        except:
            pass

        try:
            self.destroy()
        except:
            pass

        # logging.info("Forcefully terminating the process...")
        try:
            os.kill(os.getpid(), signal.SIGTERM)
        except:
            pass

    # <<< Functions to be used in top shelf of the window. >>>

    def choose_file(self): # function to choose the file from the computer
        logging.debug("choose_file button pressed!")
        self.file_path = filedialog.askopenfilename()
        self.textbox.delete(1.0, 'end')
        self.textbox.insert(1.0, self.file_path)
        logging.info(f"File selected: {self.file_path}")

    def analyze_button_top(self):

        logging.info("Analysis started!")


        if not hasattr(self, 'file_path'):
            logging.error("File path not selected!")
            messagebox.showerror("Error", "Please select a file first.")
            return

        threading.Thread(target=self.load_and_analyze).start()


        # self.load_and_analyze()

    # <<< Functions to be used in bottom shelf of the window. >>>

    def get_organism(self, selected_organism: str):
        self.selected_organism = selected_organism

        if selected_organism == 'prokaryote' or selected_organism == 'eukaryote':
            self.genome_gene_menu2.configure(values=['--select--','genome', 'gene'])
        else:
            self.genome_gene_menu2.configure(values=['--select--','genome', 'protein*'])

    def get_genome_gene(self, selected_genome_gene: str):
        self.selected_genome_gene = selected_genome_gene

        if selected_genome_gene == 'genome':
            self.accid_taxid_menu3.configure(values=['--select--','accession', 'taxon'])
        elif selected_genome_gene == 'gene':
            self.accid_taxid_menu3.configure(values=['--select--','accession', 'gene-id', 'symbol', 'taxon'])
        else:
            self.accid_taxid_menu3.configure(values=[])        

    def get_accid_taxid(self, selected_accid_taxid: str):
        self.selected_accid_taxid = selected_accid_taxid

    def download_json_button(self):
        logging.debug("download_json button pressed!")
        logging.info("Selected to download JSON file.")

        if not hasattr(self, 'selected_organism') or self.selected_organism == '--select--':
            logging.error("Organism selection is required.")
            messagebox.showerror("Error", "Please select an Organism.")
            return

        if self.selected_organism == 'eukaryote' or self.selected_organism == 'virus':
            logging.error("In this version, only prokaryotic organisms are covered.")
            messagebox.showinfo("We are sorry!", "Feature Coming Soon!")
            return

        if not hasattr(self, 'selected_genome_gene') or self.selected_genome_gene == '--select--':
            logging.error("Info About selection is required.")
            messagebox.showerror("Error", "Please select genome/gene.")
            return

        if self.selected_genome_gene == 'gene':
            logging.error("In this version, only genomes are covered.")
            messagebox.showinfo("We are sorry!", "Feature Coming Soon!")
            return
            
        if not hasattr(self, 'selected_accid_taxid') or self.selected_accid_taxid == '--select--':
            logging.error("accid_taxid selection is required.")
            messagebox.showerror("Error", "Please select a accession/taxon option.")
            return

        self.text_input_for_download = self.inputtextbox.get("1.0", "end").strip()
        print(len(self.text_input_for_download))

        if not self.text_input_for_download:
            logging.error(f"{self.selected_accid_taxid} input is required for downloading.")
            messagebox.showerror("Error!", f"Please enter the {self.selected_accid_taxid}.")
            return

        if self.selected_accid_taxid == "accession" and len(self.text_input_for_download) > 2:
            self.listed_string = self.string_to_list(self.text_input_for_download)
            self.downloaded_file_name = f"{self.listed_string[0]}" + "_and_more.json"
            logging.debug(f"setting file name to: {self.downloaded_file_name}")
        else:
            try:
                self.genus, self.species = self.text_input_for_download.strip("\"").split(" ")
                self.downloaded_file_name = f"{self.genus}" + "_" + f"{self.species}" + ".json"
                logging.debug(f"setting file name to: {self.downloaded_file_name}")
            except:
                self.single_word = self.text_input_for_download.strip("\"").strip()
                self.downloaded_file_name = f"{self.single_word}" + ".json"

        self.file_path = str(os.path.join(self.saving_file_path, self.downloaded_file_name))

        logging.info(f"Selected organism: {self.selected_organism}, Selected genome_gene: {self.selected_genome_gene}, Selected accid_taxid: {self.selected_accid_taxid}, textbox: {self.inputtextbox.get('1.0', 'end').strip()}")

        self.datasets_executable = 'datasets.exe' if os.name == 'nt' else 'datasets'

        logging.debug(f"OS name: {os.name}, datasets executable: {self.datasets_executable}")

        if self.selected_organism != 'virus':
            # self.command = f"{self.datasets_executable} summary {self.selected_genome_gene} {self.selected_accid_taxid} {self.inputtextbox.get('1.0', 'end').strip()} > ./tmp/{self.downloaded_file_name}"
            self.command = f"{self.datasets_executable} summary {self.selected_genome_gene} {self.selected_accid_taxid} {self.inputtextbox.get('1.0', 'end').strip()} > {self.file_path}"
            logging.info(f"This command will be executed to download the required JSON file: {self.command}")
        else:
            # self.command = f"{self.datasets_executable} summary {self.selected_organism} {self.selected_genome_gene} {self.selected_accid_taxid} {self.inputtextbox.get('1.0', 'end').strip()} --virus > ./tmp/{self.downloaded_file_name}"
            self.command = f"{self.datasets_executable} summary {self.selected_organism} {self.selected_genome_gene} {self.selected_accid_taxid} {self.inputtextbox.get('1.0', 'end').strip()} --virus > {self.file_path}"
            logging.info(f"This command will be executed to download the required JSON file: {self.command}")

        os.makedirs("tmp", exist_ok=True)

        
        # self.file_path = os.path.abspath(f"./tmp/{self.downloaded_file_name}")
        if os.path.exists(self.file_path):
            logging.info("File already exists!")
            messagebox.showinfo("Info", f"File already exists at {self.file_path}! You can directly go for analysis.")
            return

        if not is_internet_there():
            messagebox.showerror("No internet!", "Internet connection is required to download.")
            return

        logging.info("Internet connection available. Starting download...")

        # Start the download in a new thread
        download_thread = threading.Thread(target=self.start_downloading)
        download_thread.start()
        
        download_log_thread = threading.Thread(target=self.start_download_log)
        download_log_thread.start()
        
    def start_download_log(self):
        while True:
            # print(f"Entered the while loop. the value of self.end_download_logging is: {self.end_download_logging}")
            if not self.end_download_logging:  # Wait for the download to complete
                self.elapsed_time = time.time() - self.downloading_start_time
                logging.info(f"Downloading metadata. Time elapsed: {self.elapsed_time:.2f} seconds.")
                time.sleep(5)  # Log progress every 5 seconds
            else:
                break
    
    def start_downloading(self):
        self.end_download_logging = False
        self.downloading_start_time = time.time()
        # print(f"the value of self.end_download_logging is: {self.end_download_logging}")

        try:
            # print(f"Entered the try loop. the value of self.end_download_logging is: {self.end_download_logging}")
            self.success, self.output = pass_to_cmd(self.command)
            # print(f"the value of self.success is: {self.success}")
            if self.success:
                # print(f"Entered the success loop. the value of self.end_download_logging is: {self.end_download_logging}")
                self.end_download_logging = True
                # print(f"the value of self.end_download_logging is set to: {self.end_download_logging}")
                logging.info('JSON file downloaded successfully!')
                logging.info(f"Download complete in {self.elapsed_time:.2f} seconds.")
                if not self.file_path:
                    self.file_path = os.getcwd()
                # self.file_path = os.path.abspath('./tmp/downloaded.json')
            else:
                # print(f"Entered the else loop. the value of self.end_download_logging is: {self.end_download_logging}")
                self.end_download_logging = True
                # print(f"the value of self.end_download_logging is set to: {self.end_download_logging}")
                logging.error(f"Error downloading JSON file: {self.output}")
                messagebox.showerror("Error", f"Error downloading JSON file: {self.output}")
        except:
            # print(f"Entered the except loop. the value of self.end_download_logging is: {self.end_download_logging}")
            self.end_download_logging = True
        
        while True:
            # print(f"Entered the while loop. the value of self.end_download_logging is: {self.end_download_logging}")
            if not self.end_download_logging:  # Wait for the download to complete
                self.elapsed_time = time.time() - self.downloading_start_time
                logging.info(f"Downloading metadata. Time elapsed: {self.elapsed_time:.2f} seconds.")
                time.sleep(5)  # Log progress every 5 seconds
            else:
                break
        
    def analyze_button_bottom(self):
        logging.info("Analysis started from downloaded json file!")
        
        try:
            # self.file_path = os.path.abspath(f"./tmp/{self.downloaded_file_name}")
            # self.file_path = self.file_path
            logging.debug(f"Downloaded JSON File path: {self.file_path}")
        except:
            # self.file_path = os.path.abspath(f"./tmp/downloaded.json")
            # self.file_path = self.file_path
            logging.debug(f"Downloaded JSON File path: {self.file_path}")

        if not os.path.exists(self.file_path):
            logging.error("Downloaded JSON file not found!")
            messagebox.showerror("Error", "Please download the file first!. If already downloaded, Re-download the file and if the error still persists, rename the downloaded file to 'downloaded.json'.")
            return

        threading.Thread(target=self.load_and_analyze).start()

        # self.load_and_analyze()
    
    def string_to_list(self, string):
        if ',' in string:
            return string.split(',')
        if '\t' in string:
            return string.split('\t')
        if ' ' in string:
            return string.split(' ')
    
    def toggle_debug_mode(self):
        if self.debug_var.get() == "on":
            logging.getLogger().setLevel(logging.DEBUG)
            logging.debug("Debug mode enabled")
            messagebox.showinfo("Debug Mode", "Debug logging is now enabled.")
        else:
            logging.getLogger().setLevel(logging.INFO)
            logging.info("Debug mode disabled")
            messagebox.showinfo("Debug Mode", "Debug logging is now disabled.")

    
if __name__ == '__main__':
    # root = customtkinter.CTk()
    multiprocessing.freeze_support()
    metaminer_gui = metaminer()
    metaminer_gui.mainloop()

    metaminer_gui.cleanup()

    


