# DiagnosticDigestBio441

Meant to provide recommendations for restriction enzyme combinations when performing diagnostic digests

Input files needed:
1. .txt containing sequence of construct with gDNA insert in forward orientation
2. .csv containing enzyme names, recognition sites in one letter nomenclature, cut positions offset from end of recognition site (e.g. HindIII-HF,AAGCTT,-5,-1)
# Running the application
1. Clone repo 
2. Install Flask:
  pip install Flask
3. Run UI:
  python suggester_ui.py
 
