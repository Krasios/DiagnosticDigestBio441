from flask import Flask,request,render_template,flash
from werkzeug.utils import secure_filename
import tempfile
import os

from suggester import runSuggester


UPLOAD_DIR = tempfile.gettempdir()
ALLOWED_FILETYPES = set(['txt','csv'])

app = Flask(__name__)
app.secret_key = os.urandom(24)

def valid_file(file):
    return '.' in file and file.split('.')[-1].lower() in ALLOWED_FILETYPES
 
@app.route('/',methods=['GET','POST'])
def select_insert():
    global PARAMS
    if request.method == 'POST':
        if 'sequence' not in request.files or request.files['sequence'].filename == '':
            flash('Missing Sequence File')
            return render_template('webapp.html')

        if 'enzymes' not in request.files or request.files['enzymes'].filename == '':
            flash('Missing Enzyme File')
            return render_template('webapp.html')

        seqFile = secure_filename(request.files['sequence'].filename)
        if '.' not in seqFile or seqFile.split('.')[-1].lower() != 'txt':
            flash('Invalid Sequence File')
            return render_template('webapp.html')
        request.files['sequence'].save(os.path.join(UPLOAD_DIR,seqFile))
        enzymeFile = secure_filename(request.files['enzymes'].filename)
        if '.' not in enzymeFile or enzymeFile.split('.')[-1].lower() != 'csv':
            flash('Invalid Enzyme File')
            return render_template('webapp.html')
        request.files['enzymes'].save(os.path.join(UPLOAD_DIR,enzymeFile))
        insertPositions = (int(request.form['insertBegin']),int(request.form['insertEnd']))
        intronBeginning = [int(x) for x in request.form.getlist('intronBegin')]
        intronEnding = [int(x) for x in request.form.getlist('intronEnd')]
        intronPositions = zip(intronBeginning,intronEnding)
        circular = True if request.form['circular'] == 'True' else False
        ladder = int(request.form['ladder'])
        results = runSuggester(seqFile,enzymeFile,insertPositions,intronPositions,circular,ladder)
        return render_template('results.html', results=results)

    return render_template('webapp.html')

app.run(debug=True, port=5000)

