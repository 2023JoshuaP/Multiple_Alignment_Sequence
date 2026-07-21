import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
import time
import sys
import os
import tracemalloc

from core.fasta_parser import parse_fasta, get_sequence, get_ids
from algorithms.clustalw import ClustalW
from algorithms.muscle import MUSCLE
from algorithms.simulated_annealing import SimulatedAnnealing
from core.sp_score import sp_score

class RedirectText(object):
    def __init__(self, text_widget):
        self.output = text_widget

    def write(self, string):
        self.output.insert(tk.END, string)
        self.output.see(tk.END) # Auto-scroll al final

    def flush(self):
        pass

class MSAApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Herramienta de Alineamiento Múltiple de Secuencias (MSA)")
        self.geometry("900x700")
        
        style = ttk.Style(self)
        if "clam" in style.theme_names():
            style.theme_use("clam")
            
        self.filepath = None
        self.sequences = []
        self.ids = []
        
        self.create_widgets()
        
    def create_widgets(self):
        main_frame = ttk.Frame(self, padding="20 20 20 20")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        title_lbl = ttk.Label(main_frame, text="Múltiple Sequence Alignment", font=("Helvetica", 16, "bold"))
        title_lbl.pack(pady=(0, 20))
        
        config_frame = ttk.LabelFrame(main_frame, text="Configuración del Análisis", padding="10 10 10 10")
        config_frame.pack(fill=tk.X, pady=(0, 20))
        data_frame = ttk.Frame(config_frame)
        data_frame.pack(fill=tk.X, pady=5)
        
        self.lbl_file = ttk.Label(data_frame, text="Datos: Usando secuencias de prueba por defecto.", foreground="blue")
        self.lbl_file.pack(side=tk.LEFT, padx=(0, 10))
        
        btn_load = ttk.Button(data_frame, text="Cargar Archivo FASTA", command=self.load_fasta)
        btn_load.pack(side=tk.RIGHT)
        btn_clear = ttk.Button(data_frame, text="Usar Secuencias por Defecto", command=self.use_defaults)
        btn_clear.pack(side=tk.RIGHT, padx=5)

        algo_frame = ttk.Frame(config_frame)
        algo_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(algo_frame, text="Algoritmo a ejecutar:").pack(side=tk.LEFT, padx=(0, 15))
        
        self.algo_var = tk.StringVar(value="4")
        ttk.Radiobutton(algo_frame, text="ClustalW (Progresivo)", variable=self.algo_var, value="1").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(algo_frame, text="MUSCLE (Iterativo)", variable=self.algo_var, value="2").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(algo_frame, text="Simulated Annealing (Heurístico)", variable=self.algo_var, value="3").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(algo_frame, text="Comparar Todos", variable=self.algo_var, value="4").pack(side=tk.LEFT, padx=5)

        self.btn_run = ttk.Button(config_frame, text="Ejecutar Alineamiento", command=self.start_alignment_thread)
        self.btn_run.pack(fill=tk.X, pady=10)
        res_frame = ttk.LabelFrame(main_frame, text="Consola de Resultados", padding="10 10 10 10")
        res_frame.pack(fill=tk.BOTH, expand=True)
        
        self.text_area = tk.Text(res_frame, wrap=tk.NONE, font=("Courier", 10), bg="white", fg="black")
        vsb = ttk.Scrollbar(res_frame, orient="vertical", command=self.text_area.yview)
        hsb = ttk.Scrollbar(res_frame, orient="horizontal", command=self.text_area.xview)
        self.text_area.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        self.text_area.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        
        res_frame.grid_rowconfigure(0, weight=1)
        res_frame.grid_columnconfigure(0, weight=1)
        
        redir = RedirectText(self.text_area)
        sys.stdout = redir
        
        self.use_defaults()

    def load_fasta(self):
        filepath = filedialog.askopenfilename(
            title="Seleccionar Archivo FASTA",
            filetypes=(("FASTA files", "*.fasta *.fa"), ("Todos los archivos", "*.*"))
        )
        if filepath:
            try:
                records = parse_fasta(filepath)
                self.sequences = get_sequence(records)
                self.ids = get_ids(records)
                self.filepath = filepath
                self.lbl_file.config(text=f"Datos: {os.path.basename(filepath)} ({len(self.sequences)} secuencias)")
                print(f"Archivo cargado exitosamente: {filepath}")
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo cargar el archivo FASTA:\n{e}")

    def use_defaults(self):
        self.filepath = None
        self.sequences = [
            "MKTAYIAKQRQISFV",
            "MKTAYIAKQRQISFVKS",
            "MKTAIAKQRQISFV",
            "GGGGCCCCTTTTAAAA",
        ]
        self.ids = ["Seq1", "Seq2", "Seq3", "Seq4"]
        self.lbl_file.config(text="Datos: Usando secuencias de prueba por defecto.")
        print("Datos reseteados a las secuencias de ejemplo.")

    def format_alignment_str(self, profile: list[str], labels: list[str]) -> str:
        out = ""
        width = max(len(lbl) for lbl in labels) + 2
        for lbl, seq in zip(labels, profile):
            out += f"{lbl:<{width}}{seq}\n"
        return out

    def start_alignment_thread(self):
        self.btn_run.config(state=tk.DISABLED)
        self.text_area.delete(1.0, tk.END)
        print(f"Iniciando análisis para {len(self.sequences)} secuencias...\n")
        
        t = threading.Thread(target=self.run_alignment)
        t.daemon = True
        t.start()

    def run_alignment(self):
        opcion = self.algo_var.get()
        correr_clustal = opcion in ["1", "4"]
        correr_muscle = opcion in ["2", "4"]
        correr_sa = opcion in ["3", "4"]
        
        resultados = []

        try:
            if correr_clustal:
                print("="*50)
                print(" EJECUTANDO CLUSTALW (Progresivo - Árbol NJ) ")
                print("="*50)
                clustal = ClustalW()
                tracemalloc.start()
                start = time.time()
                msa_c, ids_c, tree_c, _ = clustal.align(self.sequences, self.ids)
                t_c = time.time() - start
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                peak_mb = peak / 1048576.0
                score_c = sp_score(msa_c)
                resultados.append(("ClustalW", t_c, peak_mb, score_c, msa_c, ids_c))
                
                os.makedirs("benchmarks", exist_ok=True)
                filepath_c = os.path.join("benchmarks", "arbol_clustalw.png")
                try:
                    from core.tree_plotter import plot_phylogenetic_tree
                    plot_phylogenetic_tree(tree_c, filepath_c, title="Árbol Guía (Neighbor-Joining)", is_nj=True)
                    print(f"-> Gráfica NJ guardada en: {filepath_c}")
                    os.system(f"xdg-open '{filepath_c}' &")
                except Exception as ex:
                    print(f"-> No se pudo graficar el árbol: {ex}")
                
            if correr_muscle:
                print("\n" + "="*50)
                print(" EJECUTANDO MUSCLE (Iterativo - Refinamiento) ")
                print("="*50)
                muscle = MUSCLE(max_iters=3)
                tracemalloc.start()
                start = time.time()
                msa_m, ids_m, score_m, tree_m = muscle.align(self.sequences, self.ids)
                t_m = time.time() - start
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                peak_mb = peak / 1048576.0
                resultados.append(("MUSCLE", t_m, peak_mb, score_m, msa_m, ids_m))

                os.makedirs("benchmarks", exist_ok=True)
                filepath_m = os.path.join("benchmarks", "arbol_muscle.png")
                try:
                    from core.tree_plotter import plot_phylogenetic_tree
                    plot_phylogenetic_tree(tree_m, filepath_m, title="Árbol Guía (UPGMA)", is_nj=False)
                    print(f"-> Gráfica UPGMA guardada en: {filepath_m}")
                    os.system(f"xdg-open '{filepath_m}' &")
                except Exception as ex:
                    print(f"-> No se pudo graficar el árbol: {ex}")

            if correr_sa:
                print("\n" + "="*50)
                print(" EJECUTANDO SIMULATED ANNEALING (Heurístico) ")
                print("="*50)
                sa = SimulatedAnnealing(initial_temp=50.0, cooling_rate=0.9, max_iters=500)
                tracemalloc.start()
                start = time.time()
                msa_s, ids_s, score_s, tree_s = sa.align(self.sequences, self.ids)
                t_s = time.time() - start
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                peak_mb = peak / 1048576.0
                resultados.append(("Simulated Annealing", t_s, peak_mb, score_s, msa_s, ids_s))
                
                if tree_s:
                    os.makedirs("benchmarks", exist_ok=True)
                    filepath_s = os.path.join("benchmarks", "arbol_sa.png")
                    try:
                        from core.tree_plotter import plot_phylogenetic_tree
                        plot_phylogenetic_tree(tree_s, filepath_s, title="Árbol Guía (Simulated Annealing)", is_nj=True)
                        print(f"-> Gráfica SA guardada en: {filepath_s}")
                        os.system(f"xdg-open '{filepath_s}' &")
                    except Exception as ex:
                        print(f"-> No se pudo graficar el árbol: {ex}")

            print("\n\n" + "#"*50)
            print(" COMPARATIVA FINAL DE RESULTADOS ")
            print("#"*50)
            
            with open("results.txt", "w", encoding="utf-8") as f:
                for nombre, tiempo, ram, score, msa, msa_ids in resultados:
                    # Imprimir en consola
                    print(f"\n[ ALGORITMO: {nombre} ]")
                    print(f"  Tiempo de ejecución : {tiempo:.4f} segundos")
                    print(f"  Uso Pico de Memoria : {ram:.2f} MB")
                    print(f"  SP-Score Objetivo   : {score:.1f} (mayor es mejor)")
                    print("-" * 50)
                    print(self.format_alignment_str(msa, msa_ids))
                    print("-" * 50)
                    
                    # Guardar en archivo log
                    f.write(f"\n[ ALGORITMO: {nombre} ]\n")
                    f.write(f"  Tiempo de ejecución : {tiempo:.4f} segundos\n")
                    f.write(f"  Uso Pico de Memoria : {ram:.2f} MB\n")
                    f.write(f"  SP-Score Objetivo   : {score:.1f} (mayor es mejor)\n")
                    f.write("-" * 50 + "\n")
                    f.write(self.format_alignment_str(msa, msa_ids) + "\n")
                    f.write("-" * 50 + "\n")
                    
            print("\n¡Análisis completado exitosamente! (Resultados guardados en results.txt)")
            
        except Exception as e:
            print(f"\nERROR DURANTE LA EJECUCIÓN:\n{e}")
        finally:
            self.after(0, lambda: self.btn_run.config(state=tk.NORMAL))

if __name__ == "__main__":
    app = MSAApp()
    app.mainloop()
