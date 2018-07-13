from flask import Flask, request, jsonify, Response
from flask_frozen import Freezer
app = Flask(__name__)
freezer = Freezer(app)

import json
import random

def make_seq(length):
    return ''.join([random.choice(list('ACTG')) for x in range(length)])

def make_feat(refseq):
    return {
        "type": "gene", "start": 5975, "end": 9744, "score": 0.84, "strand": 1,
        "name": refseq + ".au9.g1002",
        "uniqueID": refseq + ".au9.g1002",
        "subfeatures": [
            {
                "type": "mRNA", "start": 5975, "end": 9744, "score": 0.84, "strand": 1,
                "name": refseq + ".au9.g1002.t1",
                "uniqueID": refseq + ".au9.g1002.t2",
                "subfeatures": [
                    { "type": "five_prime_UTR", "start": 5975, "end": 6109, "score": 0.98, "strand": 1 },
                    { "type": "start_codon", "start": 6110, "end": 6112, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 7142, "end": 7319, "score": 1, "strand": 1, "phase": 2 },
                    { "type": "CDS",         "start": 7411, "end": 7687, "score": 1, "strand": 1, "phase": 1 },
                    { "type": "CDS",         "start": 7748, "end": 7850, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 7953, "end": 8098, "score": 1, "strand": 1, "phase": 2 },
                    { "type": "CDS",         "start": 8166, "end": 8320, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 6110, "end": 6148, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 6615, "end": 6683, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 6758, "end": 7040, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 8419, "end": 8614, "score": 1, "strand": 1, "phase": 1 },
                    { "type": "CDS",         "start": 8708, "end": 8811, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 8927, "end": 9239, "score": 1, "strand": 1, "phase": 1 },
                    { "type": "CDS",         "start": 9414, "end": 9494, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "stop_codon",  "start": 9492, "end": 9494,             "strand": 1, "phase": 0 },
                    { "type": "three_prime_UTR", "start": 9495, "end": 9744, "score": 0.86, "strand": 1 }
                ]
            },
            {
                "type": "mRNA", "start": 5975, "end": 9744, "score": 0.84, "strand": 1,
                "name": refseq + ".au9.g1002.t2",
                "uniqueID": refseq + ".au9.g1002.t2",
                "subfeatures": [
                    { "type": "five_prime_UTR", "start": 5975, "end": 6109, "score": 0.98, "strand": 1 },
                    { "type": "start_codon", "start": 6110, "end": 6112, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 6110, "end": 6148, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 6615, "end": 6683, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 6758, "end": 7040, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 8166, "end": 8320, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 8419, "end": 8614, "score": 1, "strand": 1, "phase": 1 },
                    { "type": "CDS",         "start": 8708, "end": 8811, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "CDS",         "start": 8927, "end": 9239, "score": 1, "strand": 1, "phase": 1 },
                    { "type": "CDS",         "start": 9414, "end": 9494, "score": 1, "strand": 1, "phase": 0 },
                    { "type": "stop_codon",  "start": 9492, "end": 9494,             "strand": 1, "phase": 0 },
                    { "type": "three_prime_UTR", "start": 9495, "end": 9744, "score": 0.86, "strand": 1 }
                ]
            }
        ]
    }

def feat2searchLoc(key, f):
    return {
        'name': f['name'],
        'location': {
            'ref': key,
            'start': f['start'],
            'end': f['end'],
            'tracks': [
                # ???
            ],
            'objectName': f['name'],
        }
    }

CHROM_NAMES = ('chrA', 'chrB', 'chrX')

seqs = {
    key: make_seq(random.randint(10000, 30000))
    for key in CHROM_NAMES
}

features = {
    key: [make_feat(key)]
    for key in CHROM_NAMES
}

@app.route('/refSeqs.json')
def ref_seqs():
    data = [
        {
            "length": len(seqs[key]),
            "name": key,
            "start": 0,
            "end": len(seqs[key]),
            "seqChunkSize":20000
        }
        for key in seqs.keys()
    ]
    resp = Response(
        response=json.dumps(data, indent=2),
        status=200,
        mimetype="application/json"
    )
    return resp

@app.route('/stats/global')
def stats_global():
    data = {
        # Setting this high will return "too much data to show"
        "featureDensity": 0.02,
    }
    return jsonify(**data)

@app.route('/names')
def search():
    data = []
    if request.args.get('equals'):
        for key in CHROM_NAMES:
            for f in features[key]:
                if f['name'] == request.args.get('equals'):
                    data.append(feat2searchLoc(key, f))

    elif request.args.get('startswith'):
        for key in CHROM_NAMES:
            for f in features[key]:
                if f['name'].startswith(request.args.get('startswith')):
                    data.append(feat2searchLoc(key, f))

    resp = Response(
        response=json.dumps(data),
        status=200,
        mimetype="application/json"
    )
    return resp

@app.route('/features/<refseq>')
def feats(refseq):
    if request.args.get('sequence') == 'true':
        start = int(request.args.get('start'))
        end = int(request.args.get('end'))

        data = {
            'features':[
                {
                    "seq": ''.join(seqs[refseq][start:end]),
                    "start": start,
                    "end": end
                },
            ]
        }

    else:
        data = {
            "features": features[refseq]
        }
    return jsonify(**data)

# @freezer.register_generator
# def feats():
    # for chrom in CHROM_NAMES:
        # # yield 'feats', {'refseq': chrom}
        # for i in range(0, len(seqs[chrom]), 5000):
            # yield '/features/%s?sequence=true&start=%s&end=%s' % (chrom, i, i + 5000)

@app.after_request
def apply_caching(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Origin, X-Requested-With, Content-Type, Accept"
    return response

if __name__ == "__main__":
    freezer.freeze()
    app.run(debug=True)

