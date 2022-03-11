from web_app import db


class Isomirnas(db.Model):
    isomirna = db.Column(db.Text, primary_key=True)
    ts_full_mirdb_short_comparasion = db.Column(db.Float)
    ts_short_mirdb_short_comparasion = db.Column(db.Float)

class Genes(db.Model):
    gene = db.Column(db.Text, primary_key=True)

class Targets(db.Model):
    isomirna = db.Column(db.Text, db.ForeignKey('isomirnas.isomirna', ondelete='CASCADE'), primary_key=True, index=True)
    gene = db.Column(db.Text, db.ForeignKey('genes.gene', ondelete='CASCADE'), primary_key=True, index=True)
    target_score = db.Column(db.SmallInteger)
    cwcs = db.Column(db.Float)

class Comparasion_ts_1(db.Model):
    isomirna = db.Column(db.Text, db.ForeignKey('isomirnas.isomirna', ondelete='CASCADE'), primary_key=True)
    ts_count = db.Column(db.SmallInteger)
    ts_count_1 = db.Column(db.SmallInteger)
    ts_intersection_1 = db.Column(db.Float)
    ts_union_1 = db.Column(db.Float)

class Comparasion_mirdb_1(db.Model):
    isomirna = db.Column(db.Text, db.ForeignKey('isomirnas.isomirna', ondelete='CASCADE'), primary_key=True)
    mirdb_count = db.Column(db.SmallInteger)
    mirdb_count_1 = db.Column(db.SmallInteger)
    mirdb_intersection_1 = db.Column(db.Float)
    mirdb_union_1 = db.Column(db.Float)
