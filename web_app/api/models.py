from web_app import db


class Molecules(db.Model):
    # Molecule = gene or isomiR
    molecule = db.Column(db.Text, primary_key=True)


class Cancers(db.Model):
    cancer = db.Column(db.Text, primary_key=True)


class Expression(db.Model):
    molecule = db.Column(db.Text, db.ForeignKey('molecules.molecule', ondelete='CASCADE'), primary_key=True, index=True)
    cancer = db.Column(db.Text, db.ForeignKey('cancers.cancer', ondelete='CASCADE'), primary_key=True, index=True)
    tmm = db.Column(db.ARRAY(db.Numeric(precision=4, scale=1)), nullable=False)
    tpm = db.Column(db.ARRAY(db.Numeric(precision=4, scale=1)), nullable=False)
    highly_expressed = db.Column(db.Boolean, nullable=False)
