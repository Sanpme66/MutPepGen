import tkinter as tk
import tkinter.messagebox as messagebox
import customtkinter as ctk
import os
from tkinter import filedialog
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageTk
import threading
import time
import json
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import sys

# Set appearance mode and color theme
ctk.set_appearance_mode("System")  # Default to system theme
ctk.set_default_color_theme("blue")  # Professional blue theme

# Default database location
DEFAULT_DB_PATH = os.path.join(os.getcwd(), "database")
H=1280
W=720
__version__="v1.0.0-dev"

class MutationPeptideApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        # Configure window
        self.title("Canimmune-Mutpep | Mutation-Derived Peptide Generator")
        self.geometry(f"{H}x{W}")
        # Set icon after window creation
        # self.iconphoto(False, tk.PhotoImage(file="icons/dna.png"))
        
        # self.after(100, self.set_app_icon)
        # self.set_app_icon()
        self.iconbitmap("assets/icons/dna.ico") if os.path.exists("assets/icons/dna.ico") else None
        
        # Initialize variables
        self.initialize_variables()
        
        # Configure grid layout
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure((0, 1, 2), weight=1)
        
        # Create sidebar
        self.create_sidebar()
        
        # Create main content area with tabs
        self.create_main_content()
        
        # Initialize the UI
        self.initialize_ui()
        
        # Initialize database if it doesn't exist
        self.setup_database()
        
    def set_app_icon(self):
        """Set the application icon after the window is created"""
        # For Windows and macOS compatibility
        try:
            if os.name == 'nt':  # Windows
                if os.path.exists("assets/icons/canimmune_icon.ico"):
                    self.iconbitmap("assets/icons/canimmune_icon.ico")
            else:  # macOS or Linux
                if os.path.exists("assets/icons/canimmune_icon.png"):
                    icon_img = tk.PhotoImage(file="assets/icons/canimmune_icon.png")
                    self.iconphoto(True, icon_img)
        except Exception as e:
            print(f"Could not set application icon: {e}")
            
    def initialize_variables(self):
        """Initialize application variables"""
        self.input_files = []
        self.current_file = None
        self.column_mapping = {"enst_id": None, "mutation": None}
        self.available_columns = []
        self.has_selected_columns = False
        self.database_path = os.path.join(os.getcwd(), "database", "ensembl_sequences.fasta")
        self.output_dir = os.path.join(os.getcwd(), "results")
        self.peptide_window = tk.IntVar(value=25)
        self.include_sequence_info = tk.BooleanVar(value=True)
        self.num_threads = tk.IntVar(value=4)
        self.version = __version__
        self.processing_in_progress = False
        self.df = None
        self.sequence_db = {}
        
        # Set up color scheme
        self.colors = {
            "primary": "#1976D2",
            "primary_dark": "#0D47A1", 
            "secondary": "#00BFA5",
            "secondary_dark": "#00796B",
            "accent": "#FF9800",
            "accent_dark": "#EF6C00",
            "text_light": "#F5F5F5",
            "text_dark": "#212121",
            "bg_light": "#FAFAFA",
            "bg_dark": "#263238"
        }

    def create_sidebar(self):
        """Create sidebar with application controls"""
        # Main sidebar frame
        self.sidebar_frame = ctk.CTkFrame(self, width=280, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(20, weight=1)
        
        # App logo and title
        logo_frame = ctk.CTkFrame(self.sidebar_frame, fg_color="transparent")
        logo_frame.grid(row=0, column=0, padx=20, pady=(20, 10), sticky="ew")
        
        try:
            # Load and display logo if available
            self.logo_image = ctk.CTkImage(
                light_image=Image.open("icons/canimmune_logo.png"),
                dark_image=Image.open("icons/canimmune_logo.png"),
                size=(40, 40)
            )
            logo_label = ctk.CTkLabel(logo_frame, image=self.logo_image, text="")
            logo_label.pack(side="left", padx=(0, 10))
        except:
            pass
            
        title_label = ctk.CTkLabel(
            logo_frame, 
            text="Canimmune-MutPep", 
            font=ctk.CTkFont(size=24, weight="bold")
        )
        title_label.pack(side="left")
        
        subtitle_label = ctk.CTkLabel(
            self.sidebar_frame,
            text="Mutation-Derived Peptide Generator",
            font=ctk.CTkFont(size=12)
        )
        subtitle_label.grid(row=1, column=0, padx=20, pady=(0, 15))
        
        # Separator
        separator1 = ctk.CTkFrame(self.sidebar_frame, height=1)
        separator1.grid(row=2, column=0, sticky="ew", padx=15, pady=(0, 15))
        
        # Section label - Input Data
        input_section_label = ctk.CTkLabel(
            self.sidebar_frame,
            text="INPUT DATA",
            font=ctk.CTkFont(size=12, weight="bold")
        )
        input_section_label.grid(row=3, column=0, padx=20, pady=(0, 10), sticky="w")
        
        # Input file buttons
        self.file_button = ctk.CTkButton(
            self.sidebar_frame, 
            text="Select Mutation File", 
            command=self.select_input_file,
            height=36
        )
        self.file_button.grid(row=4, column=0, padx=20, pady=(0, 10), sticky="ew")
        
        self.column_button = ctk.CTkButton(
            self.sidebar_frame, 
            text="Map Columns", 
            command=self.show_column_mapping, 
            state="disabled",
            height=36
        )
        self.column_button.grid(row=5, column=0, padx=20, pady=(0, 10), sticky="ew")
        
        # Section label - Database
        db_section_label = ctk.CTkLabel(
            self.sidebar_frame,
            text="SEQUENCE DATABASE",
            font=ctk.CTkFont(size=12, weight="bold")
        )
        db_section_label.grid(row=6, column=0, padx=20, pady=(15, 10), sticky="w")
        
        self.db_button = ctk.CTkButton(
            self.sidebar_frame, 
            text="Select Sequence Database", 
            command=self.select_database,
            height=36
        )
        self.db_button.grid(row=7, column=0, padx=20, pady=(0, 10), sticky="ew")
        
        # Section label - Parameters
        param_section_label = ctk.CTkLabel(
            self.sidebar_frame,
            text="ANALYSIS PARAMETERS",
            font=ctk.CTkFont(size=12, weight="bold")
        )
        param_section_label.grid(row=8, column=0, padx=20, pady=(15, 10), sticky="w")
        
        # Peptide window size frame
        peptide_frame = ctk.CTkFrame(self.sidebar_frame, fg_color="transparent")
        peptide_frame.grid(row=9, column=0, padx=20, pady=(0, 10), sticky="ew")
        
        peptide_label = ctk.CTkLabel(
            peptide_frame, 
            text="Peptide Window Size:"
        )
        peptide_label.grid(row=0, column=0, padx=0, pady=(0, 5), sticky="w")
        
        self.peptide_value = ctk.CTkLabel(
            peptide_frame, 
            textvariable=self.peptide_window
        )
        self.peptide_value.grid(row=0, column=1, padx=5, pady=(0, 5), sticky="e")
        
        self.peptide_slider = ctk.CTkSlider(
            peptide_frame, 
            from_=5, 
            to=21, 
            number_of_steps=8, 
            variable=self.peptide_window,
            command=self.update_peptide_window
        )
        self.peptide_slider.grid(row=1, column=0, columnspan=2, pady=(0, 5), sticky="ew")
        
        # Output options frame
        output_options_frame = ctk.CTkFrame(self.sidebar_frame, fg_color="transparent")
        output_options_frame.grid(row=10, column=0, padx=20, pady=(0, 10), sticky="ew")
        
        self.seq_info_checkbox = ctk.CTkCheckBox(
            output_options_frame,
            text="Include sequence info in headers",
            variable=self.include_sequence_info
        )
        self.seq_info_checkbox.grid(row=0, column=0, sticky="w")
        
        # Threads frame
        threads_frame = ctk.CTkFrame(self.sidebar_frame, fg_color="transparent")
        threads_frame.grid(row=11, column=0, padx=20, pady=(0, 10), sticky="ew")
        
        threads_label = ctk.CTkLabel(
            threads_frame, 
            text="Processing Threads:"
        )
        threads_label.grid(row=0, column=0, padx=0, pady=(0, 5), sticky="w")
        
        self.threads_value = ctk.CTkLabel(
            threads_frame, 
            textvariable=self.num_threads
        )
        self.threads_value.grid(row=0, column=1, padx=5, pady=(0, 5), sticky="e")
        
        self.threads_slider = ctk.CTkSlider(
            threads_frame, 
            from_=1, 
            to=16, 
            number_of_steps=15, 
            variable=self.num_threads
        )
        self.threads_slider.grid(row=1, column=0, columnspan=2, pady=(0, 5), sticky="ew")
        
        # Output directory
        self.output_button = ctk.CTkButton(
            self.sidebar_frame, 
            text="Set Results Directory", 
            command=self.select_output_dir,
            height=36
        )
        self.output_button.grid(row=12, column=0, padx=20, pady=(0, 10), sticky="ew")
        
        # Separator
        separator2 = ctk.CTkFrame(self.sidebar_frame, height=1)
        separator2.grid(row=13, column=0, sticky="ew", padx=15, pady=(15, 15))
        
        # Run button at bottom
        self.run_button = ctk.CTkButton(
            self.sidebar_frame, 
            text="Generate Peptides", 
            command=self.run_analysis, 
            height=48,
            fg_color=self.colors["secondary"],
            hover_color=self.colors["secondary_dark"],
            font=ctk.CTkFont(size=16, weight="bold")
        )
        self.run_button.grid(row=14, column=0, padx=20, pady=(0, 20), sticky="ew")
        
        # Status indicators
        self.status_frame = ctk.CTkFrame(self.sidebar_frame)
        self.status_frame.grid(row=15, column=0, padx=20, pady=(0, 15), sticky="ew")
        
        self.status_label = ctk.CTkLabel(
            self.status_frame,
            text="Ready",
            font=ctk.CTkFont(size=12)
        )
        self.status_label.pack(padx=10, pady=10)
        
        # Version info
        self.version_label = ctk.CTkLabel(
            self.sidebar_frame,
            text=f"MutPep {self.version} \n © Chen Li Lab / Purcell Lab \n Monash University",
            font=ctk.CTkFont(size=10)
        )
        self.version_label.grid(row=16, column=0, padx=20, pady=(10, 5))
        
        # Help and appearance settings
        bottom_frame = ctk.CTkFrame(self.sidebar_frame, fg_color="transparent")
        bottom_frame.grid(row=17, column=0, padx=20, pady=(0, 15), sticky="ew")
        
        self.help_button = ctk.CTkButton(
            bottom_frame, 
            text="Help", 
            command=self.show_help,
            width=80,
            height=28
        )
        self.help_button.grid(row=0, column=0, padx=(0, 5), sticky="w")
        
        appearance_label = ctk.CTkLabel(
            bottom_frame, 
            text="Theme:", 
            font=ctk.CTkFont(size=11)
        )
        appearance_label.grid(row=0, column=1, padx=(10, 5), sticky="e")
        
        self.appearance_menu = ctk.CTkOptionMenu(
            bottom_frame, 
            values=["Light", "Dark", "System"],
            command=self.change_appearance_mode_event,
            width=80,
            height=28
        )
        self.appearance_menu.grid(row=0, column=2, sticky="e")
        self.appearance_menu.set("System")

    def create_main_content(self):
        """Create main content area with tabs"""
        # Create tabview for main content
        self.tabview = ctk.CTkTabview(self)
        self.tabview.grid(row=0, column=1, rowspan=4, padx=(20, 20), pady=(20, 20), sticky="nsew")
        
        # Create tabs
        self.tab_dashboard = self.tabview.add("Dashboard")
        self.tab_data = self.tabview.add("Data Explorer")
        self.tab_log = self.tabview.add("Processing Log")
        self.tab_results = self.tabview.add("Results")
        
        # Configure each tab's grid
        for tab in [self.tab_dashboard, self.tab_data, self.tab_log, self.tab_results]:
            tab.grid_columnconfigure(0, weight=1)
            tab.grid_rowconfigure(1, weight=1)
        
        # Setup dashboard tab content
        self.setup_dashboard_tab()
        
        # Setup data explorer tab
        self.setup_data_explorer_tab()
        
        # Setup log tab content
        self.setup_log_tab()
        
        # Setup results tab content
        self.setup_results_tab()

    def setup_dashboard_tab(self):
        """Setup dashboard tab content"""
        # Dashboard header
        dashboard_title = ctk.CTkLabel(
            self.tab_dashboard,
            text="Mutation-Derived Peptide Database Generator",
            font=ctk.CTkFont(size=22, weight="bold")
        )
        dashboard_title.grid(row=0, column=0, padx=20, pady=(20, 5), sticky="w")
        
        # adding Workflow png
        # Load the image
        try:
            image_path = "assets/icons/workflow_fig.png"
            image = ctk.CTkImage(light_image=Image.open(image_path), size=(int(H*0.7), (W*0.3)))  # Adjust size as needed

            # Create a label with the image
            dashboard_image = ctk.CTkLabel(self.tab_dashboard, image=image, text="")
            dashboard_image.grid(row=0, column=0, padx=20, pady=(20, 10), sticky="w")
        except Exception as e:
            self.log_message(f"Error loading Workflow image: {str(e)}", "error")

        # Subtitle below the image
        dashboard_subtitle = ctk.CTkLabel(
            self.tab_dashboard,
            text="Generate peptide sequences centered on mutation sites for LC-MS/MS database search",
            font=ctk.CTkFont(size=14)
        )
        dashboard_subtitle.grid(row=1, column=0, padx=20, pady=(0, 20), sticky="w")
        
        
        # Dashboard content frame
        dashboard_frame = ctk.CTkFrame(self.tab_dashboard)
        dashboard_frame.grid(row=2, column=0, padx=20, pady=(0, 20), sticky="nsew")
        dashboard_frame.grid_columnconfigure(0, weight=1)
        dashboard_frame.grid_columnconfigure(1, weight=1)
        
        # Status section - Left Column
        self.status_left = ctk.CTkFrame(dashboard_frame)
        self.status_left.grid(row=0, column=0, padx=(20, 10), pady=20, sticky="nsew")
        self.status_left.grid_columnconfigure(0, weight=1)
        
        status_title = ctk.CTkLabel(
            self.status_left,
            text="Analysis Status",
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color=self.colors["primary"],
            text_color=self.colors["text_light"],
            corner_radius=6
        )
        status_title.grid(row=0, column=0, padx=15, pady=(15, 15), sticky="ew")
        
        # Create status indicators
        self.input_file_status = self.create_status_indicator(self.status_left, "Mutation Data File", 1)
        self.column_status = self.create_status_indicator(self.status_left, "Column Mapping", 2)
        self.database_status = self.create_status_indicator(self.status_left, "Sequence Database", 3)
        self.peptide_status = self.create_status_indicator(self.status_left, "Peptide Window Size", 4)
        self.output_status = self.create_status_indicator(self.status_left, "Results Directory", 5)
        
        # Dashboard right column - Quick stats and visualization
        self.status_right = ctk.CTkFrame(dashboard_frame)
        self.status_right.grid(row=0, column=1, padx=(10, 20), pady=20, sticky="nsew")
        self.status_right.grid_columnconfigure(0, weight=1)
        
        workflow_title = ctk.CTkLabel(
            self.status_right,
            text="Workflow Overview",
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color=self.colors["primary"],
            text_color=self.colors["text_light"],
            corner_radius=6
        )
        workflow_title.grid(row=0, column=0, padx=15, pady=(15, 15), sticky="ew")
        
        # Create a workflow diagram
        self.workflow_frame = ctk.CTkFrame(self.status_right)
        self.workflow_frame.grid(row=1, column=0, padx=15, pady=10, sticky="nsew")
        
        workflow_steps = [
            {"step": "1. Load Mutation Data", "description": "Upload CSV, TSV, or MAF file with mutation information"},
            {"step": "2. Map Data Columns", "description": "Identify ENST_ID and mutation information columns"}, 
            {"step": "3. Configure Parameters", "description": "Set peptide window size and output options"},
            {"step": "4. Process Mutations", "description": "Extract sequences and generate peptides for each mutation"},
            {"step": "5. Export Results", "description": "Save peptide sequences as FASTA and analysis summary"}
        ]
        
        for i, step_info in enumerate(workflow_steps):
            step_frame = ctk.CTkFrame(self.workflow_frame)
            step_frame.grid(row=i, column=0, padx=10, pady=5, sticky="ew")
            
            step_label = ctk.CTkLabel(
                step_frame,
                text=step_info["step"],
                font=ctk.CTkFont(size=14, weight="bold")
            )
            step_label.grid(row=0, column=0, padx=10, pady=(10, 5), sticky="w")
            
            desc_label = ctk.CTkLabel(
                step_frame,
                text=step_info["description"],
                font=ctk.CTkFont(size=12)
            )
            desc_label.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="w")

    def setup_data_explorer_tab(self):
        """Setup data explorer tab content"""
        # Data explorer header
        data_title = ctk.CTkLabel(
            self.tab_data,
            text="Data Explorer",
            font=ctk.CTkFont(size=20, weight="bold")
        )
        data_title.grid(row=0, column=0, padx=20, pady=(20, 5), sticky="w")
        
        data_subtitle = ctk.CTkLabel(
            self.tab_data,
            text="Explore and validate your mutation data before processing",
            font=ctk.CTkFont(size=14)
        )
        data_subtitle.grid(row=1, column=0, padx=20, pady=(0, 20), sticky="w")
        
        # Data explorer main frame
        data_frame = ctk.CTkFrame(self.tab_data)
        data_frame.grid(row=2, column=0, padx=20, pady=(0, 20), sticky="nsew")
        data_frame.grid_columnconfigure(0, weight=1)
        data_frame.grid_rowconfigure(1, weight=1)
        
        # Controls frame
        controls_frame = ctk.CTkFrame(data_frame, fg_color="transparent")
        controls_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        
        refresh_button = ctk.CTkButton(
            controls_frame,
            text="Refresh Data",
            command=self.refresh_data_view,
            width=120
        )
        refresh_button.pack(side="left", padx=(10, 5))
        
        self.data_info_label = ctk.CTkLabel(
            controls_frame,
            text="No data loaded",
            font=ctk.CTkFont(size=12)
        )
        self.data_info_label.pack(side="left", padx=10)
        
        # Data table frame (with scrollbars)
        table_frame = ctk.CTkFrame(data_frame)
        table_frame.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="nsew")
        table_frame.grid_rowconfigure(0, weight=1)
        table_frame.grid_columnconfigure(0, weight=1)
        
        # Create scrollable frame for the data table
        self.data_view = ctk.CTkTextbox(table_frame, wrap="none")
        self.data_view.grid(row=0, column=0, sticky="nsew")
        
        # Add horizontal scrollbar
        h_scrollbar = ctk.CTkScrollbar(data_frame, orientation="horizontal", command=self.data_view.xview)
        h_scrollbar.grid(row=2, column=0, sticky="ew", padx=10)
        self.data_view.configure(xscrollcommand=h_scrollbar.set)

    def setup_log_tab(self):
        """Setup log tab content"""
        # Log header
        log_title = ctk.CTkLabel(
            self.tab_log,
            text="Processing Log",
            font=ctk.CTkFont(size=20, weight="bold")
        )
        log_title.grid(row=0, column=0, padx=20, pady=(20, 5), sticky="w")
        
        # Log frame
        log_frame = ctk.CTkFrame(self.tab_log)
        log_frame.grid(row=1, column=0, padx=20, pady=(10, 20), sticky="nsew")
        log_frame.grid_rowconfigure(0, weight=1)
        log_frame.grid_columnconfigure(0, weight=1)
        
        # Log textbox with custom styling
        self.log_textbox = ctk.CTkTextbox(log_frame, width=500)
        self.log_textbox.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        # Configure tags for different message types
        # self.log_textbox.tag_config("header", font=ctk.CTkFont(size=14, weight="bold"))
        # self.log_textbox.tag_config("subheader", font=ctk.CTkFont(size=13, weight="bold"))
        # self.log_textbox.tag_config("info")
        # self.log_textbox.tag_config("warning")
        # self.log_textbox.tag_config("error")
        # self.log_textbox.tag_config("success")
        self.log_textbox.tag_config("header", foreground="#1976D2")
        self.log_textbox.tag_config("subheader", foreground="#0277BD")
        self.log_textbox.tag_config("info", foreground="#212121")
        self.log_textbox.tag_config("warning", foreground="#EF6C00")
        self.log_textbox.tag_config("error", foreground="#C62828")
        self.log_textbox.tag_config("success", foreground="#00695C")
    def setup_results_tab(self):
        """Setup results tab content"""
        # Results header
        results_title = ctk.CTkLabel(
            self.tab_results,
            text="Analysis Results",
            font=ctk.CTkFont(size=20, weight="bold")
        )
        results_title.grid(row=0, column=0, padx=20, pady=(20, 5), sticky="w")
        
        # Results main frame with two columns
        results_frame = ctk.CTkFrame(self.tab_results)
        results_frame.grid(row=1, column=0, padx=20, pady=(10, 20), sticky="nsew")
        results_frame.grid_columnconfigure(0, weight=1)
        results_frame.grid_columnconfigure(1, weight=1)
        results_frame.grid_rowconfigure(0, weight=1)
        
        # Left column - Results textbox
        left_frame = ctk.CTkFrame(results_frame)
        left_frame.grid(row=0, column=0, padx=(10, 5), pady=10, sticky="nsew")
        left_frame.grid_rowconfigure(1, weight=1)
        left_frame.grid_columnconfigure(0, weight=1)
        
        stats_label = ctk.CTkLabel(
            left_frame,
            text="Result Statistics",
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color=self.colors["primary"],
            text_color=self.colors["text_light"],
            corner_radius=6
        )
        stats_label.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        
        self.results_textbox = ctk.CTkTextbox(left_frame)
        self.results_textbox.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="nsew")
        self.results_textbox.insert("0.0", "No analysis results available yet.\n\nPlease configure the analysis parameters and click 'Generate Peptides' to start processing.")
        
        # Right column - Visualization
        right_frame = ctk.CTkFrame(results_frame)
        right_frame.grid(row=0, column=1, padx=(5, 10), pady=10, sticky="nsew")
        right_frame.grid_rowconfigure(1, weight=1)
        right_frame.grid_columnconfigure(0, weight=1)
        
        viz_label = ctk.CTkLabel(
            right_frame,
            text="Visualization",
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color=self.colors["primary"],
            text_color=self.colors["text_light"],
            corner_radius=6
        )
        viz_label.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        
        # Frame for visualization that will be populated after analysis
        self.viz_frame = ctk.CTkFrame(right_frame)
        self.viz_frame.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="nsew")
        
        # Placeholder text for visualization
        placeholder_label = ctk.CTkLabel(
            self.viz_frame,
            text="Visualizations will appear here after analysis",
            font=ctk.CTkFont(size=14)
        )
        placeholder_label.place(relx=0.5, rely=0.5, anchor=tk.CENTER)
        
        # Bottom frame for export options
        export_frame = ctk.CTkFrame(self.tab_results, fg_color="transparent")
        export_frame.grid(row=2, column=0, padx=20, pady=(0, 20), sticky="ew")
        
        export_button = ctk.CTkButton(
            export_frame,
            text="Export All Results",
            command=self.export_all_results,
            height=36,
            fg_color=self.colors["secondary"],
            hover_color=self.colors["secondary_dark"]
        )
        export_button.pack(side="left", padx=10)
        
        export_fasta_button = ctk.CTkButton(
            export_frame,
            text="Export FASTA Only",
            command=self.export_fasta_only,
            height=36
        )
        export_fasta_button.pack(side="left", padx=10)
        
        export_summary_button = ctk.CTkButton(
            export_frame,
            text="Export Summary Report",
            command=self.export_summary_report,
            height=36
        )
        export_summary_button.pack(side="left", padx=10)
        
    def initialize_ui(self):
        """Initialize the UI with welcome message and default values"""
        self.log_message("MutPep | Mutation-Derived Peptide Generator", "header")
        self.log_message("=========================================", "header")
        self.log_message("\nWelcome to MutPep!", "subheader")
        self.log_message("\nThis tool helps generate peptide sequences centered on mutation sites for:")
        self.log_message("  • Neoantigen prediction and vaccine development")
        self.log_message("  • Immunogenicity studies")
        self.log_message("  • Epitope mapping")
        self.log_message("  • Structural biology analysis")
        
        self.log_message("\nQuick Start Guide:", "subheader")
        self.log_message("  1. Upload your mutation data file (CSV, TSV, or MAF)")
        self.log_message("  2. Map the required columns (ENST_ID and mutation information)")
        self.log_message("  3. Set your desired peptide window size")
        self.log_message("  4. Click 'Generate Peptides' to process your data")
        
        # Initialize status indicators
        self.update_status(self.database_status, True, f"Default database ✓") 
        self.update_peptide_window()
        self.update_output_status()

    def setup_database(self):
        """Check if the database exists, create if not"""
        if not os.path.exists(DEFAULT_DB_PATH):
            os.makedirs(DEFAULT_DB_PATH, exist_ok=True)
            
        db_file = os.path.join(DEFAULT_DB_PATH, "ensembl_sequences.fasta")
        
        if not os.path.exists(db_file):
            self.log_message("\nNotice: Default sequence database not found.", "warning")
            self.log_message("You will need to provide a FASTA file containing ENST sequences.", "info")
            self.update_status(self.database_status, False)
        else:
            self.log_message(f"Found default sequence database at {db_file}", "info")
            self.load_sequence_database(db_file)
            
    def load_sequence_database(self, db_path):
        """Load the sequence database into memory"""
        try:
            self.log_message(f"Loading sequence database from {db_path}...", "info")
            start_time = time.time()
            
            # Clear existing database
            self.sequence_db = {}
            
            # Load sequences from FASTA file
            for record in SeqIO.parse(db_path, "fasta"):
                # Extract ENST ID from the record ID (assuming format like "ENST00000123456.1")
                enst_id = record.id.split('.')[0]
                self.sequence_db[enst_id] = str(record.seq)
            
            end_time = time.time()
            self.log_message(f"Loaded {len(self.sequence_db)} sequences in {end_time - start_time:.2f} seconds", "success")
            
            # Update database status
            db_name = os.path.basename(db_path)
            self.update_status(self.database_status, True, f"{db_name} ✓")
        except Exception as e:
            self.log_message(f"Error loading sequence database: {str(e)}", "error")
            self.update_status(self.database_status, False)
            
    def create_status_indicator(self, parent, label_text, row, required=True):
        """Create a status indicator for the dashboard"""
        frame = ctk.CTkFrame(parent, fg_color="transparent")
        frame.grid(row=row, column=0, padx=15, pady=(5, 10), sticky="ew")
        
        label = ctk.CTkLabel(
            frame, 
            text=f"{label_text}:",
            font=ctk.CTkFont(size=13)
        )
        label.grid(row=0, column=0, padx=0, pady=0, sticky="w")
        
        status_text = "Required" if required else "Optional"
        status_color = "#C62828" if required else "#78909C"
        
        status = ctk.CTkLabel(
            frame, 
            text=status_text, 
            text_color=status_color,
            font=ctk.CTkFont(size=13)
        )
        status.grid(row=0, column=1, padx=10, pady=0, sticky="e")
        
        return status
    
    def update_status(self, status_label, is_set=False, text=None, optional=False):
        """Update a status indicator with styling"""
        if is_set:
            display_text = text if text else "Set ✓"
            status_label.configure(text=display_text, text_color="#00695C")
        else:
            if optional:
                status_label.configure(text="Optional", text_color="#78909C")
            else:
                status_label.configure(text="Required", text_color="#C62828")
                
    def update_peptide_window(self, value=None):
        """Update peptide window size when slider changes"""
        window_size = self.peptide_window.get()
        self.update_status(self.peptide_status, True, f"{window_size} amino acids ✓")
        
    def change_appearance_mode_event(self, new_appearance_mode):
        """Change app appearance mode"""
        ctk.set_appearance_mode(new_appearance_mode)
        
    def log_message(self, message, tag=None):
        """Add a message to the log with optional tag styling"""
        self.log_textbox.configure(state="normal")
        if tag:
            self.log_textbox.insert("end", message + "\n", tag)
        else:
            self.log_textbox.insert("end", message + "\n")
        self.log_textbox.configure(state="normal")
        self.log_textbox.see("end")
        
    def select_input_file(self):
        """Select input mutation file"""
        filepath = filedialog.askopenfilename(
            title="Select Mutation Data File",
            filetypes=[
                ("Mutation Files", "*.csv;*.tsv;*.maf"),
                ("CSV Files", "*.csv"),
                ("TSV Files", "*.tsv"),
                ("MAF Files", "*.maf"),
                ("All Files", "*.*")
            ]
        )
        
        if filepath:
            # Check if it's a valid file type
            ext = os.path.splitext(filepath)[1].lower()
            if ext in ['.csv', '.tsv', '.maf']:
                self.input_files = [filepath]
                self.current_file = filepath
                file_name = os.path.basename(filepath)
                
                self.update_status(self.input_file_status, True, f"{file_name} ✓")
                self.log_message(f"Selected mutation file: {file_name}", "info")
                
                # Reset column mapping
                self.column_mapping = {"enst_id": None, "mutation": None}
                self.update_status(self.column_status, False)
                self.has_selected_columns = False
                
                # Enable column mapping button
                self.column_button.configure(state="normal")
                
                # Try to load the file
                self.load_file(filepath)
            else:
                self.update_status(self.input_file_status, False)
                messagebox.showwarning("Invalid File", "Please select a CSV, TSV, or MAF file")
        else:
            # User canceled selection
            if not self.input_files:
                self.update_status(self.input_file_status, False)
                
    def load_file(self, filepath):
        """Load the selected file and detect columns"""
        try:
            ext = os.path.splitext(filepath)[1].lower()
            self.log_message(f"Loading file {filepath}...", "info")
            
            if ext == '.csv':
                self.df = pd.read_csv(filepath)
            elif ext == '.tsv':
                self.df = pd.read_csv(filepath, sep='\t')
            elif ext == '.maf':
                # For MAF files, skip lines starting with #
                self.df = pd.read_csv(filepath, sep='\t', comment='#')
            
            # Update data explorer
            self.refresh_data_view()
            
            # Get column names
            self.available_columns = list(self.df.columns)
            
            self.log_message(f"Successfully loaded file with {len(self.df)} rows and {len(self.available_columns)} columns", "success")
            
        except Exception as e:
            self.log_message(f"Error loading file: {str(e)}", "error")
            self.df = None
            
    def refresh_data_view(self):
        """Refresh the data view in the Data Explorer tab"""
        self.data_view.configure(state="normal")
        self.data_view.delete("0.0", "end")
        
        if self.df is not None:
            # Update info label
            self.data_info_label.configure(text=f"Rows: {len(self.df)} | Columns: {len(self.df.columns)}")
            
            # Display first 100 rows in the data view
            max_rows = min(100, len(self.df))
            preview_df = self.df.head(max_rows).copy()
            
            # Format as a table with borders
            table_text = preview_df.to_string(index=True)
            self.data_view.insert("0.0", table_text)
            
            # Add a note if showing only a subset
            if len(self.df) > 100:
                self.data_view.insert("end", f"\n\nNote: Showing first 100 of {len(self.df)} rows")
        else:
            self.data_info_label.configure(text="No data loaded")
            self.data_view.insert("0.0", "No data loaded. Please select a file first.")
        
        self.data_view.configure(state="disabled")
        
    def show_column_mapping(self):
        """Display dialog to map columns from the file"""
        if not self.input_files or self.df is None:
            messagebox.showerror("Error", "Please select a valid mutation file first")
            return
            
        # Create a new dialog window
        dialog = ctk.CTkToplevel(self)
        dialog.title("Map Data Columns")
        dialog.geometry("550x500")
        dialog.transient(self)  # Make dialog modal
        dialog.grab_set()
        
        # Configure dialog grid
        dialog.grid_columnconfigure(0, weight=1)
        
        # Create header
        header_label = ctk.CTkLabel(
            dialog,
            text="Column Mapping for Mutation Data",
            font=ctk.CTkFont(size=18, weight="bold"),
            text_color=self.colors["primary"]
        )
        header_label.grid(row=0, column=0, padx=20, pady=(20, 5), sticky="w")
        
        instruction_label = ctk.CTkLabel(
            dialog,
            text="Select the columns that correspond to the required data fields:",
            font=ctk.CTkFont(size=13)
        )
        instruction_label.grid(row=1, column=0, padx=20, pady=(0, 15), sticky="w")
        
        # Create frame for mapping fields
        mapping_frame = ctk.CTkFrame(dialog)
        mapping_frame.grid(row=2, column=0, padx=20, pady=(0, 20), sticky="nsew")
        
        # Required field - Transcript ID (ENST)
        enst_label = ctk.CTkLabel(
            mapping_frame,
            text="Transcript ID (ENST):",
            font=ctk.CTkFont(weight="bold"),
            text_color=self.colors["primary"]
        )
        enst_label.grid(row=0, column=0, padx=15, pady=(15, 5), sticky="w")
        
        enst_var = tk.StringVar()
        enst_dropdown = ctk.CTkOptionMenu(
            mapping_frame,
            values=self.available_columns,
            variable=enst_var,
            dynamic_resizing=False,
            width=200
        )
        enst_dropdown.grid(row=0, column=1, padx=15, pady=(15, 5), sticky="e")
        
        # Try to auto-select based on common column names
        for col in self.available_columns:
            col_lower = col.lower()
            if "transcript" in col_lower or "enst" in col_lower or "feature" in col_lower:
                enst_var.set(col)
                break
        
        enst_desc = ctk.CTkLabel(
            mapping_frame,
            text="The column containing Ensembl transcript IDs\n(e.g., ENST00000000000)",
            font=ctk.CTkFont(size=11)
        )
        enst_desc.grid(row=1, column=0, columnspan=2, padx=15, pady=(0, 15), sticky="w")
        
        # Required field - Mutation
        mutation_label = ctk.CTkLabel(
            mapping_frame,
            text="Mutation Description:",
            font=ctk.CTkFont(weight="bold"),
            text_color=self.colors["primary"]
        )
        mutation_label.grid(row=2, column=0, padx=15, pady=(0, 5), sticky="w")
        
        mutation_var = tk.StringVar()
        mutation_dropdown = ctk.CTkOptionMenu(
            mapping_frame,
            values=self.available_columns,
            variable=mutation_var,
            dynamic_resizing=False,
            width=200
        )
        mutation_dropdown.grid(row=2, column=1, padx=15, pady=(0, 5), sticky="e")
        
        # Try to auto-select based on common column names
        for col in self.available_columns:
            col_lower = col.lower()
            if "hgvs" in col_lower or "mutation" in col_lower or "variant" in col_lower or "amino_acid" in col_lower:
                mutation_var.set(col)
                break
        
        mutation_desc = ctk.CTkLabel(
            mapping_frame,
            text="The column containing mutation information\n(e.g., p.V600E or c.1799T>A)",
            font=ctk.CTkFont(size=11)
        )
        mutation_desc.grid(row=3, column=0, columnspan=2, padx=15, pady=(0, 15), sticky="w")
        
        # Preview section
        preview_label = ctk.CTkLabel(
            dialog,
            text="Data Preview (First 5 Rows):",
            font=ctk.CTkFont(size=14, weight="bold"),
            text_color=self.colors["primary"]
        )
        preview_label.grid(row=3, column=0, padx=20, pady=(0, 5), sticky="w")
        
        preview_frame = ctk.CTkFrame(dialog)
        preview_frame.grid(row=4, column=0, padx=20, pady=(0, 15), sticky="nsew")
        
        preview_text = ctk.CTkTextbox(preview_frame, height=130, width=500)
        preview_text.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        # Load preview data
        if self.df is not None:
            preview_df = self.df.head(5)
            preview_text.insert("0.0", preview_df.to_string(index=False))
        else:
            preview_text.insert("0.0", "No data available for preview")
        
        # Buttons
        button_frame = ctk.CTkFrame(dialog, fg_color="transparent")
        button_frame.grid(row=5, column=0, padx=20, pady=(0, 20), sticky="e")
        
        def save_mapping():
            # Save the column mapping
            if not enst_var.get() or not mutation_var.get():
                messagebox.showerror("Error", "Please select both required fields")
                return
                
            self.column_mapping["enst_id"] = enst_var.get()
            self.column_mapping["mutation"] = mutation_var.get()
            self.has_selected_columns = True
            
            # Update UI status
            mapping_summary = f"{enst_var.get()} & {mutation_var.get()} ✓"
            self.update_status(self.column_status, True, mapping_summary)
            self.log_message(f"Column mapping set - Transcript ID: {enst_var.get()}, Mutation: {mutation_var.get()}", "info")
            
            dialog.destroy()
        
        def cancel():
            dialog.destroy()
        
        save_button = ctk.CTkButton(
            button_frame,
            text="Save Mapping",
            command=save_mapping,
            fg_color=self.colors["secondary"],
            hover_color=self.colors["secondary_dark"]
        )
        save_button.grid(row=0, column=0, padx=(0, 10))
        
        cancel_button = ctk.CTkButton(
            button_frame,
            text="Cancel",
            command=cancel
        )
        cancel_button.grid(row=0, column=1)
        
    def select_database(self):
        """Select reference protein database file"""
        db_path = filedialog.askopenfilename(
            title="Select Sequence Database",
            filetypes=[
                ("FASTA Files", "*.fasta;*.fa"),
                ("All Files", "*.*")
            ]
        )
        
        if db_path:
            self.database_path = db_path
            self.load_sequence_database(db_path)
        else:
            # User canceled, keep current database
            pass
            
    def select_output_dir(self):
        """Select output directory for results"""
        output_dir = filedialog.askdirectory(
            title="Select Output Directory"
        )
        
        if output_dir:
            self.output_dir = output_dir
            self.update_output_status()
            self.log_message(f"Results will be saved to: {output_dir}", "info")
            
    def update_output_status(self):
        """Update output directory status"""
        dir_name = os.path.basename(self.output_dir) or "results"
        if not os.path.exists(self.output_dir):
            try:
                os.makedirs(self.output_dir)
                self.log_message(f"Created output directory: {self.output_dir}", "info")
            except Exception as e:
                self.log_message(f"Error creating output directory: {str(e)}", "error")
                
        self.update_status(self.output_status, True, f"{dir_name} ✓")
        
    def validate_inputs(self):
        """Validate all required inputs before running analysis"""
        missing = []
        
        if not self.input_files or self.df is None:
            missing.append("Input mutation file")
        
        if not self.has_selected_columns:
            missing.append("Column mapping")
            
        if not self.sequence_db:
            missing.append("Sequence database")
        
        if missing:
            error_message = "Missing required inputs:\n" + "\n".join([f"- {item}" for item in missing])
            messagebox.showerror("Missing Inputs", error_message)
            return False
            
        return True
        
    def run_analysis(self):
        """Run the peptide generation analysis"""
        if self.processing_in_progress:
            messagebox.showinfo("Processing", "Analysis is already running. Please wait.")
            return
            
        if not self.validate_inputs():
            return
            
        # Set processing flag
        self.processing_in_progress = True
        self.status_label.configure(text="Processing...")
        
        # Switch to the log tab to show progress
        self.tabview.set("Processing Log")
        
        # Clear previous results
        self.results_textbox.delete("0.0", "end")
        
        # Log analysis start
        self.log_message("\n" + "="*50, "header")
        self.log_message("Starting Mutation Peptide Generation", "header")
        self.log_message("="*50, "header")
        self.log_message(f"Time: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.log_message(f"Input File: {os.path.basename(self.current_file)}")
        self.log_message(f"Peptide Window Size: {self.peptide_window.get()} amino acids")
        self.log_message(f"Include sequence info in headers: {'Yes' if self.include_sequence_info.get() else 'No'}")
        self.log_message(f"Processing threads: {self.num_threads.get()}")
        self.log_message(f"Results Directory: {self.output_dir}")
        self.log_message("-"*50)
        
        # Start processing in a separate thread
        threading.Thread(target=self.process_mutations).start()
        
    def process_mutations(self):
        """Process mutations and generate peptides"""
        try:
            window_size = self.peptide_window.get()
            half_window = window_size // 2
            
            self.log_message("Step 1: Preparing mutation data...", "subheader")
            
            # Get the columns we need
            enst_column = self.column_mapping["enst_id"]
            mutation_column = self.column_mapping["mutation"]
            
            # Create results directory if it doesn't exist
            os.makedirs(self.output_dir, exist_ok=True)
            
            # Initialize counters
            total_mutations = len(self.df)
            processed_mutations = 0
            successful_peptides = 0
            failed_peptides = 0
            
            # Prepare results dictionary
            results = {
                "mutation_peptides": [],
                "stats": {
                    "total_mutations": total_mutations,
                    "processed_mutations": 0,
                    "successful_peptides": 0,
                    "failed_peptides": 0,
                    "invalid_transcripts": 0,
                    "invalid_mutations": 0
                }
            }
            
            # Create output FASTA file
            fasta_path = os.path.join(self.output_dir, "mutation_peptides.fasta")
            summary_path = os.path.join(self.output_dir, "analysis_summary.json")
            
            self.log_message(f"Processing {total_mutations} mutations...", "info")
            
            with open(fasta_path, 'w') as fasta_out:
                # Process each mutation
                for index, row in self.df.iterrows():
                    try:
                        # Get transcript ID
                        transcript_id = str(row[enst_column]).strip()
                        
                        # Ensure ENST format
                        if not transcript_id.startswith("ENST"):
                            transcript_id = f"ENST{transcript_id}" if transcript_id.isdigit() else transcript_id
                            
                        # Remove version number if present
                        if "." in transcript_id:
                            transcript_id = transcript_id.split(".")[0]
                        
                        # Get mutation information
                        mutation_info = str(row[mutation_column]).strip()
                        
                        # Log progress every 100 mutations
                        if index % 100 == 0 or index == total_mutations - 1:
                            self.log_message(f"Processing mutation {index+1}/{total_mutations}: {transcript_id} {mutation_info}", "info")
                        
                        # Check if transcript exists in database
                        if transcript_id not in self.sequence_db:
                            self.log_message(f"Warning: Transcript {transcript_id} not found in database", "warning")
                            results["stats"]["invalid_transcripts"] += 1
                            failed_peptides += 1
                            continue
                        
                        # Get the sequence
                        sequence = self.sequence_db[transcript_id]
                        
                        # Parse the mutation
                        try:
                            # Handling different mutation formats
                            if mutation_info.startswith("p."):
                                # Protein mutation format (e.g., p.V600E)
                                mutation_info = mutation_info[2:]  # Remove p. prefix
                                
                                # Extract position and mutation
                                position = ""
                                for char in mutation_info[1:]:
                                    if char.isdigit():
                                        position += char
                                    else:
                                        break
                                        
                                if not position:
                                    raise ValueError(f"Could not extract position from {mutation_info}")
                                    
                                position = int(position) - 1  # Convert to 0-based index
                                mutant_aa = mutation_info[-1]
                                
                                # Validate position
                                if position < 0 or position >= len(sequence):
                                    raise ValueError(f"Position {position+1} is out of range for sequence length {len(sequence)}")
                                
                                # Extract peptide region
                                start = max(0, position - half_window)
                                end = min(len(sequence), position + half_window + 1)
                                
                                peptide = sequence[start:end]
                                
                                # Create mutant peptide by replacing the amino acid at mutation position
                                rel_pos = position - start
                                if 0 <= rel_pos < len(peptide):
                                    mutant_peptide = peptide[:rel_pos] + mutant_aa + peptide[rel_pos+1:]
                                else:
                                    raise ValueError(f"Relative position {rel_pos} is out of range for peptide length {len(peptide)}")
                                
                                # Create FASTA header
                                if self.include_sequence_info.get():
                                    header = f">{transcript_id}|{mutation_info}|pos:{position+1}|window:{window_size}|mutant"
                                else:
                                    header = f">{transcript_id}_{mutation_info}_mutant"
                                
                                # Write to FASTA file
                                fasta_out.write(f"{header}\n{mutant_peptide}\n")
                                
                                # Store in results
                                results["mutation_peptides"].append({
                                    "transcript_id": transcript_id,
                                    "mutation": mutation_info,
                                    "position": position + 1,
                                    "peptide": mutant_peptide,
                                    "original_aa": sequence[position],
                                    "mutant_aa": mutant_aa
                                })
                                
                                successful_peptides += 1
                                results["stats"]["successful_peptides"] += 1
                                
                            else:
                                # Other mutation formats not handled yet
                                self.log_message(f"Unrecognized mutation format: {mutation_info}", "warning")
                                results["stats"]["invalid_mutations"] += 1
                                failed_peptides += 1
                                continue
                                
                        except Exception as e:
                            self.log_message(f"Error processing mutation {mutation_info}: {str(e)}", "error")
                            results["stats"]["invalid_mutations"] += 1
                            failed_peptides += 1
                            continue
                            
                        processed_mutations += 1
                        results["stats"]["processed_mutations"] += 1
                        
                    except Exception as e:
                        self.log_message(f"Error processing row {index}: {str(e)}", "error")
                        failed_peptides += 1
                        continue
            
            # Update final statistics
            results["stats"]["processed_mutations"] = processed_mutations
            results["stats"]["successful_peptides"] = successful_peptides
            results["stats"]["failed_peptides"] = failed_peptides
            
            # Save summary to JSON
            with open(summary_path, 'w') as json_out:
                json.dump(results, json_out, indent=2)
            
            # Log completion
            self.log_message("\nAnalysis completed!", "header")
            self.log_message(f"Generated {successful_peptides} peptides from {processed_mutations} mutations", "success")
            self.log_message(f"Failed to process {failed_peptides} mutations", "info")
            self.log_message(f"Results saved to: {self.output_dir}", "info")
            
            # Update results tab
            self.display_results(results)
            
        except Exception as e:
            self.log_message(f"Error during analysis: {str(e)}", "error")
        finally:
            self.processing_in_progress = False
            self.status_label.configure(text="Ready")
    
    def display_results(self, results):
        """Display results in the results tab"""
        # Switch to results tab
        self.tabview.set("Results")
        
        # Update results textbox
        self.results_textbox.delete("0.0", "end")
        
        # Format results summary
        summary = f"""
# Mutation Peptide Analysis Results

## Summary Statistics
- **Total Mutations**: {results['stats']['total_mutations']}
- **Successfully Processed**: {results['stats']['processed_mutations']}
- **Peptides Generated**: {results['stats']['successful_peptides']}
- **Failed Mutations**: {results['stats']['failed_peptides']}
- **Invalid Transcripts**: {results['stats']['invalid_transcripts']}
- **Invalid Mutation Format**: {results['stats']['invalid_mutations']}

## Parameters Used
- **Peptide Window Size**: {self.peptide_window.get()} amino acids
- **Include Sequence Info**: {'Yes' if self.include_sequence_info.get() else 'No'}

## File Locations
- **FASTA File**: {os.path.join(self.output_dir, 'mutation_peptides.fasta')}
- **Analysis Summary**: {os.path.join(self.output_dir, 'analysis_summary.json')}
        """
        
        self.results_textbox.insert("0.0", summary)
        
        # Create visualization
        self.create_results_visualization(results)
    
    def create_results_visualization(self, results):
        """Create visualizations for the results"""
        # Clear existing widgets in viz_frame
        for widget in self.viz_frame.winfo_children():
            widget.destroy()
        
        # Create figure with two subplots
        fig = plt.Figure(figsize=(8, 6), dpi=100)
        
        # Add a subplot for mutation processing results
        ax1 = fig.add_subplot(121)
        labels = ['Successful', 'Failed - Invalid Transcript', 'Failed - Invalid Mutation']
        sizes = [
            results['stats']['successful_peptides'],
            results['stats']['invalid_transcripts'],
            results['stats']['invalid_mutations']
        ]
        colors = ['#2196F3', '#FFC107', '#F44336']
        ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax1.set_title('Mutation Processing Results')
        
        # Add a subplot for peptide length distribution if we have successful peptides
        if results['mutation_peptides']:
            ax2 = fig.add_subplot(122)
            peptide_lengths = [len(peptide['peptide']) for peptide in results['mutation_peptides']]
            ax2.hist(peptide_lengths, bins=10, color='#4CAF50')
            ax2.set_xlabel('Peptide Length')
            ax2.set_ylabel('Count')
            ax2.set_title('Peptide Length Distribution')
        
        # Tight layout
        fig.tight_layout()
        
        # Embed the figure in the viz_frame
        canvas = FigureCanvasTkAgg(fig, master=self.viz_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def export_all_results(self):
        """Export all results to the output directory"""
        if not os.path.exists(os.path.join(self.output_dir, "mutation_peptides.fasta")):
            messagebox.showinfo("No Results", "Please run the analysis first to generate results")
            return
            
        # All files should already be in the output directory
        messagebox.showinfo("Export Complete", f"All results have been saved to:\n{self.output_dir}")
        
        # Open the output directory in file explorer
        try:
            if os.name == 'nt':  # Windows
                os.startfile(self.output_dir)
            elif os.name == 'posix':  # macOS or Linux
                if sys.platform == 'darwin':  # macOS
                    subprocess.call(['open', self.output_dir])
                else:  # Linux
                    subprocess.call(['xdg-open', self.output_dir])
        except:
            pass
    
    def export_fasta_only(self):
        """Export only the FASTA file to a user-selected location"""
        fasta_path = os.path.join(self.output_dir, "mutation_peptides.fasta")
        if not os.path.exists(fasta_path):
            messagebox.showinfo("No Results", "Please run the analysis first to generate results")
            return
            
        # Ask for a new file location
        new_path = filedialog.asksaveasfilename(
            title="Save FASTA File",
            defaultextension=".fasta",
            filetypes=[("FASTA Files", "*.fasta"), ("All Files", "*.*")],
            initialfile="mutation_peptides.fasta"
        )
        
        if new_path:
            # Copy the file
            try:
                with open(fasta_path, 'r') as src, open(new_path, 'w') as dst:
                    dst.write(src.read())
                messagebox.showinfo("Export Complete", f"FASTA file exported to:\n{new_path}")
            except Exception as e:
                messagebox.showerror("Export Error", f"Error exporting FASTA file: {str(e)}")
    
    def export_summary_report(self):
        """Generate and export a detailed summary report"""
        json_path = os.path.join(self.output_dir, "analysis_summary.json")
        if not os.path.exists(json_path):
            messagebox.showinfo("No Results", "Please run the analysis first to generate results")
            return
            
        # Ask for a new file location
        new_path = filedialog.asksaveasfilename(
            title="Save Summary Report",
            defaultextension=".html",
            filetypes=[("HTML Files", "*.html"), ("All Files", "*.*")],
            initialfile="mutation_analysis_report.html"
        )
        
        if new_path:
            try:
                # Load the results
                with open(json_path, 'r') as f:
                    results = json.load(f)
                
                # Generate an HTML report
                html_content = self.generate_html_report_v2(results)
                
                # Write the HTML report
                with open(new_path, 'w') as f:
                    f.write(html_content)
                    
                messagebox.showinfo("Export Complete", f"Summary report exported to:\n{new_path}")
                
                # Try to open the HTML file in the default browser
                try:
                    import webbrowser
                    webbrowser.open(new_path)
                except:
                    pass
                    
            except Exception as e:
                messagebox.showerror("Export Error", f"Error generating summary report: {str(e)}")
    
    def generate_html_report(self, results):
        """Generate an HTML report from the results"""
        # Get the current date and time
        now = time.strftime("%Y-%m-%d %H:%M:%S")
        
        # Start building the HTML content
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Mutation Peptide Analysis Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}
        h1, h2, h3 {{
            color: #1976D2;
        }}
        .container {{
            background-color: #f9f9f9;
            border-radius: 8px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #1976D2;
            color: white;
        }}
        tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        .stats {{
            display: flex;
            justify-content: space-between;
            flex-wrap: wrap;
        }}
        .stat-box {{
            background-color: #fff;
            border-radius: 8px;
            padding: 15px;
            margin: 10px;
            flex: 1;
            min-width: 200px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            text-align: center;
        }}
        .stat-box h3 {{
            margin: 0;
            font-size: 16px;
        }}
        .stat-box .value {{
            font-size: 28px;
            font-weight: bold;
            color: #1976D2;
            margin: 10px 0;
        }}
        .footer {{
            margin-top: 40px;
            font-size: 12px;
            color: #777;
            text-align: center;
        }}
    </style>
</head>
<body>
    <h1>Mutation Peptide Analysis Report</h1>
    <p><strong>Generated on:</strong> {now}</p>
    
    <div class="container">
        <h2>Summary Statistics</h2>
        <div class="stats">
            <div class="stat-box">
                <h3>Total Mutations</h3>
                <div class="value">{results['stats']['total_mutations']}</div>
            </div>
            <div class="stat-box">
                <h3>Processed Mutations</h3>
                <div class="value">{results['stats']['processed_mutations']}</div>
            </div>
            <div class="stat-box">
                <h3>Peptides Generated</h3>
                <div class="value">{results['stats']['successful_peptides']}</div>
            </div>
            <div class="stat-box">
                <h3>Failed Mutations</h3>
                <div class="value">{results['stats']['failed_peptides']}</div>
            </div>
        </div>
        
        <h3>Failure Details</h3>
        <div class="stats">
            <div class="stat-box">
                <h3>Invalid Transcripts</h3>
                <div class="value">{results['stats']['invalid_transcripts']}</div>
            </div>
            <div class="stat-box">
                <h3>Invalid Mutation Format</h3>
                <div class="value">{results['stats']['invalid_mutations']}</div>
            </div>
        </div>
    </div>
    
    <div class="container">
        <h2>Analysis Parameters</h2>
        <ul>
            <li><strong>Peptide Window Size:</strong> {self.peptide_window.get()} amino acids</li>
            <li><strong>Include Sequence Info:</strong> {'Yes' if self.include_sequence_info.get() else 'No'}</li>
            <li><strong>Processing Threads:</strong> {self.num_threads.get()}</li>
            <li><strong>Input File:</strong> {os.path.basename(self.current_file)}</li>
        </ul>
    </div>
"""
        
        # Add peptide table if there are peptides
        if results['mutation_peptides']:
            # Get a sample of peptides (maximum 100 for performance)
            sample_size = min(100, len(results['mutation_peptides']))
            sample_peptides = results['mutation_peptides'][:sample_size]
            
            html += f"""
    <div class="container">
        <h2>Generated Peptides</h2>
        <p>Showing {sample_size} out of {len(results['mutation_peptides'])} peptides</p>
        <table>
            <thead>
                <tr>
                    <th>#</th>
                    <th>Transcript ID</th>
                    <th>Mutation</th>
                    <th>Position</th>
                    <th>Original AA</th>
                    <th>Mutant AA</th>
                    <th>Peptide</th>
                </tr>
            </thead>
            <tbody>
"""
            
            for i, peptide in enumerate(sample_peptides):
                html += f"""
                <tr>
                    <td>{i+1}</td>
                    <td>{peptide['transcript_id']}</td>
                    <td>{peptide['mutation']}</td>
                    <td>{peptide['position']}</td>
                    <td>{peptide['original_aa']}</td>
                    <td>{peptide['mutant_aa']}</td>
                    <td style="font-family: monospace;">{peptide['peptide']}</td>
                </tr>
"""
            
            html += """
            </tbody>
        </table>
    </div>
"""
        
        # Close the HTML
        html += f"""
    <div class="footer">
        <p>Generated by MutPep {self.version} &copy; {time.strftime("%Y")}</p>
    </div>
</body>
</html>
"""
        
        return html
    
    def generate_html_report_v2(self, results):
        """Generate an HTML report from the results"""
        # Get the current date and time
        now = time.strftime("%Y-%m-%d %H:%M:%S")
        
        # Define colors for status indicators
        success_color = "#00695C"  # Success green
        warning_color = "#FF9800"  # Warning orange
        error_color = "#C62828"    # Error red
        
        # Start building the HTML content
        html = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>CANIMMUNE-MutPep | database construction Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
            }}
            h1, h2, h3 {{
                color: #2B7DE1;
            }}
            .container {{
                background-color: #f9f9f9;
                border-radius: 8px;
                padding: 20px;
                margin-bottom: 20px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            .header {{
                display: flex;
                align-items: center;
                margin-bottom: 20px;
            }}
            .header img {{
                height: 60px;
                margin-right: 20px;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 12px;
                text-align: left;
            }}
            th {{
                background-color: #2B7DE1;
                color: white;
            }}
            tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}
            .stats {{
                display: flex;
                justify-content: space-between;
                flex-wrap: wrap;
            }}
            .stat-box {{
                background-color: #fff;
                border-radius: 8px;
                padding: 15px;
                margin: 10px;
                flex: 1;
                min-width: 200px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                text-align: center;
            }}
            .stat-box h3 {{
                margin: 0;
                font-size: 16px;
            }}
            .stat-box .value {{
                font-size: 28px;
                font-weight: bold;
                color: #2B7DE1;
                margin: 10px 0;
            }}
            .status-success {{
                color: {success_color};
            }}
            .status-warning {{
                color: {warning_color};
            }}
            .status-error {{
                color: {error_color};
            }}
            .footer {{
                margin-top: 40px;
                font-size: 12px;
                color: #777;
                text-align: center;
                border-top: 1px solid #eee;
                padding-top: 20px;
            }}
            .logo {{
                text-align: center;
                margin-bottom: 20px;
            }}
        </style>
    </head>
    <body>
        <div class="logo">
            <img src="assets/icons/canimmune_logo.png" alt="CanImmune Logo" height="80">
        </div>
        <h1>CANIMMUNE-MutPep: Cancer Neoantigen Mutation database construction Report</h1>
        <p><strong>Generated on:</strong> {now}</p>
        
        <div class="container">
            <h2>Summary Statistics</h2>
            <div class="stats">
                <div class="stat-box">
                    <h3>Total Mutations</h3>
                    <div class="value">{results['stats']['total_mutations']}</div>
                </div>
                <div class="stat-box">
                    <h3>Processed Mutations</h3>
                    <div class="value" class="status-success">{results['stats']['processed_mutations']}</div>
                </div>
                <div class="stat-box">
                    <h3>Peptides Generated</h3>
                    <div class="value" class="status-success">{results['stats']['successful_peptides']}</div>
                </div>
                <div class="stat-box">
                    <h3>Failed Mutations</h3>
                    <div class="value" class="{
                        'status-success' if results['stats']['failed_peptides'] == 0 else 
                        'status-warning' if results['stats']['failed_peptides'] < results['stats']['total_mutations'] * 0.2 else 
                        'status-error'
                    }">{results['stats']['failed_peptides']}</div>
                </div>
            </div>
            
            <h3>Failure Details</h3>
            <div class="stats">
                <div class="stat-box">
                    <h3>Invalid Transcripts</h3>
                    <div class="value" class="{
                        'status-success' if results['stats']['invalid_transcripts'] == 0 else 
                        'status-warning' if results['stats']['invalid_transcripts'] < results['stats']['total_mutations'] * 0.1 else 
                        'status-error'
                    }">{results['stats']['invalid_transcripts']}</div>
                </div>
                <div class="stat-box">
                    <h3>Invalid Mutation Format</h3>
                    <div class="value" class="{
                        'status-success' if results['stats']['invalid_mutations'] == 0 else 
                        'status-warning' if results['stats']['invalid_mutations'] < results['stats']['total_mutations'] * 0.1 else 
                        'status-error'
                    }">{results['stats']['invalid_mutations']}</div>
                </div>
            </div>
        </div>
        
        <div class="container">
            <h2>Analysis Parameters</h2>
            <ul>
                <li><strong>Peptide Window Size:</strong> {self.peptide_window.get()} amino acids</li>
                <li><strong>Include Sequence Info:</strong> {'Yes' if self.include_sequence_info.get() else 'No'}</li>
                <li><strong>Processing Threads:</strong> {self.num_threads.get()}</li>
                <li><strong>Input File:</strong> {os.path.basename(self.current_file)}</li>
            </ul>
        </div>
    """
            
        # Add peptide table if there are peptides
        if results['mutation_peptides']:
            # Get a sample of peptides (maximum 100 for performance)
            sample_size = min(100, len(results['mutation_peptides']))
            sample_peptides = results['mutation_peptides'][:sample_size]
            
            html += f"""
    <div class="container">
        <h2>Generated Peptides</h2>
        <p>Showing {sample_size} out of {len(results['mutation_peptides'])} peptides</p>
        <table>
            <thead>
                <tr>
                    <th>#</th>
                    <th>Transcript ID</th>
                    <th>Mutation</th>
                    <th>Position</th>
                    <th>Original AA</th>
                    <th>Mutant AA</th>
                    <th>Peptide</th>
                </tr>
            </thead>
            <tbody>
    """
            
            for i, peptide in enumerate(sample_peptides):
                html += f"""
                <tr>
                    <td>{i+1}</td>
                    <td>{peptide['transcript_id']}</td>
                    <td>{peptide['mutation']}</td>
                    <td>{peptide['position']}</td>
                    <td>{peptide['original_aa']}</td>
                    <td>{peptide['mutant_aa']}</td>
                    <td style="font-family: monospace;">{peptide['peptide']}</td>
                </tr>
    """
            
            html += """
            </tbody>
        </table>
    </div>
    """
        
        # Close the HTML with updated footer
        html += f"""
    <div class="footer">
        <p>Generated by CanImmune {self.version} &copy; {time.strftime("%Y")}</p>
        <p>Developed at Chen Li Lab / Purcell Lab - Monash University</p>
    </div>
    </body>
    </html>
    """
        
        return html
    
    def show_help(self):
        """Show help information"""
        help_dialog = ctk.CTkToplevel(self)
        help_dialog.title("MutPep Help")
        help_dialog.geometry("650x550")
        help_dialog.transient(self)
        help_dialog.grab_set()
        
        # Configure dialog grid
        help_dialog.grid_columnconfigure(0, weight=1)
        help_dialog.grid_rowconfigure(2, weight=1)
        
        # Create header
        header_label = ctk.CTkLabel(
            help_dialog,
            text="MutPep Help Guide",
            font=ctk.CTkFont(size=20, weight="bold"),
            text_color=self.colors["primary"]
        )
        header_label.grid(row=0, column=0, padx=20, pady=(20, 5), sticky="w")
        
        version_label = ctk.CTkLabel(
            help_dialog,
            text=f"Version {self.version}",
            font=ctk.CTkFont(size=12)
        )
        version_label.grid(row=1, column=0, padx=20, pady=(0, 15), sticky="w")
        
        # Help content in a scrollable frame
        help_frame = ctk.CTkScrollableFrame(help_dialog)
        help_frame.grid(row=2, column=0, padx=20, pady=(0, 20), sticky="nsew")
        help_frame.grid_columnconfigure(0, weight=1)
        
        # Introduction section
        intro_label = ctk.CTkLabel(
            help_frame,
            text="Introduction",
            font=ctk.CTkFont(size=16, weight="bold"),
            text_color=self.colors["primary"]
        )
        intro_label.grid(row=0, column=0, padx=10, pady=(10, 5), sticky="w")
        
        intro_text = ctk.CTkLabel(
            help_frame,
            text="MutPep is a tool for generating peptide sequences centered around "
                 "mutation sites in proteins. These peptides can be used for "
                 "neoantigen prediction, vaccine development, immunogenicity studies, "
                 "and other applications in cancer research and immunotherapy.",
            font=ctk.CTkFont(size=12),
            wraplength=580,
            justify="left"
        )
        intro_text.grid(row=1, column=0, padx=10, pady=(0, 15), sticky="w")
        
        # Workflow section
        workflow_label = ctk.CTkLabel(
            help_frame,
            text="Workflow",
            font=ctk.CTkFont(size=16, weight="bold"),
            text_color=self.colors["primary"]
        )
        workflow_label.grid(row=2, column=0, padx=10, pady=(10, 5), sticky="w")
        
        workflow_steps = [
            "1. Load mutation data (CSV, TSV, or MAF format)",
            "2. Map columns to required data fields (ENST ID and mutation)",
            "3. Select or provide a sequence database containing protein sequences",
            "4. Configure peptide window size and other parameters",
            "5. Run analysis to generate peptides for each mutation",
            "6. Export and visualize results"
        ]
        
        workflow_text = ctk.CTkLabel(
            help_frame,
            text="\n".join(workflow_steps),
            font=ctk.CTkFont(size=12),
            wraplength=580,
            justify="left"
        )
        workflow_text.grid(row=3, column=0, padx=10, pady=(0, 15), sticky="w")
        
        # Input Files section
        input_label = ctk.CTkLabel(
            help_frame,
            text="Input Files",
            font=ctk.CTkFont(size=16, weight="bold"),
            text_color=self.colors["primary"]
        )
        input_label.grid(row=4, column=0, padx=10, pady=(10, 5), sticky="w")
        
        input_text = ctk.CTkLabel(
            help_frame,
            text="MutPep accepts the following file formats:\n"
                 "• CSV/TSV: Tabular data with mutation information\n"
                 "• MAF: Mutation Annotation Format files\n\n"
                 "Required data fields:\n"
                 "• Transcript ID (ENST): Ensembl transcript identifier\n"
                 "• Mutation: Mutation description in HGVS format (e.g., p.V600E)",
            font=ctk.CTkFont(size=12),
            wraplength=580,
            justify="left"
        )
        input_text.grid(row=5, column=0, padx=10, pady=(0, 15), sticky="w")
        
        # Sequence Database section
        db_label = ctk.CTkLabel(
            help_frame,
            text="Sequence Database",
            font=ctk.CTkFont(size=16, weight="bold"),
            text_color=self.colors["primary"]
        )
        db_label.grid(row=6, column=0, padx=10, pady=(10, 5), sticky="w")
        
        db_text = ctk.CTkLabel(
            help_frame,
            text="MutPep requires a FASTA file containing protein sequences "
                 "with Ensembl transcript IDs (ENST) as headers. The application "
                 "includes a default database, but you can also provide your own.\n\n"
                 "The sequence database is used to extract the protein sequence "
                 "corresponding to each transcript ID in the mutation data.",
            font=ctk.CTkFont(size=12),
            wraplength=580,
            justify="left"
        )
        db_text.grid(row=7, column=0, padx=10, pady=(0, 15), sticky="w")
        
        # Parameters section
        params_label = ctk.CTkLabel(
            help_frame,
            text="Analysis Parameters",
            font=ctk.CTkFont(size=16, weight="bold"),
            text_color=self.colors["primary"]
        )
        params_label.grid(row=8, column=0, padx=10, pady=(10, 5), sticky="w")
        
        params_text = ctk.CTkLabel(
            help_frame,
            text="• Peptide Window Size: Length of peptides to generate "
                 "(centered on the mutation site)\n"
                 "• Include Sequence Info: Add additional information to FASTA "
                 "headers (position, window size)\n"
                 "• Processing Threads: Number of parallel threads to use for analysis",
            font=ctk.CTkFont(size=12),
            wraplength=580,
            justify="left"
        )
        params_text.grid(row=9, column=0, padx=10, pady=(0, 15), sticky="w")
        
        # Results section
        results_label = ctk.CTkLabel(
            help_frame,
            text="Results",
            font=ctk.CTkFont(size=16, weight="bold"),
            text_color=self.colors["primary"]
        )
        results_label.grid(row=10, column=0, padx=10, pady=(10, 5), sticky="w")
        
        results_text = ctk.CTkLabel(
            help_frame,
            text="Results are provided in multiple formats:\n"
                 "• Summary statistics in the application\n"
                 "• FASTA file with mutant peptide sequences\n"
                 "• JSON file with detailed analysis data\n"
                 "• HTML report with visualizations and peptide tables\n\n"
                 "You can export these results to a location of your choice "
                 "for further analysis.",
            font=ctk.CTkFont(size=12),
            wraplength=580,
            justify="left"
        )
        results_text.grid(row=11, column=0, padx=10, pady=(0, 15), sticky="w")
        
        # Support section
        support_label = ctk.CTkLabel(
            help_frame,
            text="Support",
            font=ctk.CTkFont(size=16, weight="bold"),
            text_color=self.colors["primary"]
        )
        support_label.grid(row=12, column=0, padx=10, pady=(10, 5), sticky="w")
        
        support_text = ctk.CTkLabel(
            help_frame,
            text="For additional help or to report issues, please contact:\n"
                 "chen.li@monash.edu\n\n"
                 "Documentation and tutorials are available at:\n"
                 "https://mutpep-tools.github.io/docs",
            font=ctk.CTkFont(size=12),
            wraplength=580,
            justify="left"
        )
        support_text.grid(row=13, column=0, padx=10, pady=(0, 15), sticky="w")
        
        # Close button
        close_button = ctk.CTkButton(
            help_dialog,
            text="Close",
            command=help_dialog.destroy,
            width=100
        )
        close_button.grid(row=3, column=0, padx=20, pady=(0, 20))
    
def show_splash_screen():
    """Show a splash screen while the application loads"""
    splash_root = tk.Tk()
    splash_root.overrideredirect(True)
    
    # Center the splash screen
    width = 500
    height = 300
    screen_width = splash_root.winfo_screenwidth()
    screen_height = splash_root.winfo_screenheight()
    x = (screen_width - width) // 2
    y = (screen_height - height) // 2
    splash_root.geometry(f"{width}x{height}+{x}+{y}")
    
    # Try to load the logo
    try:
        logo_path = "assets/icons/canimmune_logo.png"
        if os.path.exists(logo_path):
            logo_img = Image.open(logo_path)
            logo_img = logo_img.resize((300, 100), Image.LANCZOS)
            logo_photo = ImageTk.PhotoImage(logo_img)
            
            logo_label = tk.Label(splash_root, image=logo_photo)
            logo_label.pack(pady=20)
    except:
        # Fallback to text if image loading fails
        logo_label = tk.Label(splash_root, text="CanImmune", font=("Arial", 24, "bold"))
        logo_label.pack(pady=20)
    
    # Add a loading message
    loading_label = tk.Label(splash_root, text="Loading...", font=("Arial", 12))
    loading_label.pack(pady=10)
    
    # Add a version and copyright notice
    version_label = tk.Label(splash_root, text=f"Version 2.0.0", font=("Arial", 10))
    version_label.pack(pady=5)
    
    copyright_label = tk.Label(splash_root, text="© Chen Li Lab / Purcell Lab - Monash University", font=("Arial", 8))
    copyright_label.pack(pady=5)
    
    splash_root.update()
    
    return splash_root
def setup_application():
    """
    Prepare the application directories and resources
    """
    # Create necessary directories
    os.makedirs(DEFAULT_DB_PATH, exist_ok=True)
    os.makedirs("assets/icons", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    
    # Copy sample data files if they don't exist
    sample_db_path = os.path.join(DEFAULT_DB_PATH, "sample_sequences.fasta")
    if not os.path.exists(sample_db_path) and os.path.exists("resources/sample_sequences.fasta"):
        import shutil
        shutil.copy("resources/sample_sequences.fasta", sample_db_path)
    
    # Copy sample mutations file if it doesn't exist
    sample_mutations_path = os.path.join(os.getcwd(), "sample_mutations.csv")
    if not os.path.exists(sample_mutations_path) and os.path.exists("resources/sample_mutations.csv"):
        import shutil
        shutil.copy("resources/sample_mutations.csv", sample_mutations_path)
    
    # Check if the icons directory contains necessary files
    if not os.path.exists("assets/icons/canimmune_logo.png"):
        print("Warning: CanImmune logo not found. The application will use text headers instead.")

if __name__ == "__main__":

    # Show splash screen
    # splash = show_splash_screen()
    
    # Setup application resources
    # setup_application()
    
    # Destroy splash after a brief delay
    # splash.after(2000, splash.destroy)
    # time.sleep(2)
    app = MutationPeptideApp()
    app.mainloop()