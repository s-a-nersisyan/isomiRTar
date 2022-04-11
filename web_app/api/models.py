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


class Targets_raw(db.Model):
    isomir = db.Column(db.Text, db.ForeignKey('molecules.molecule', ondelete='CASCADE'), primary_key=True, index=True)
    target = db.Column(db.Text, db.ForeignKey('molecules.molecule', ondelete='CASCADE'), primary_key=True, index=True)
    cancer = db.Column(db.Text, db.ForeignKey('cancers.cancer', ondelete='CASCADE'), primary_key=True, index=True)
    
    mirdb_score = db.Column(db.SmallInteger)
    targetscan_score = db.Column(db.Numeric(precision=3, scale=1))
    
    spearman_corr = db.Column(db.Numeric(precision=3, scale=2))
    spearman_p_value = db.Column(db.Float)
    isomir_median_tpm = db.Column(db.Numeric(precision=4, scale=1), nullable=False)
    target_median_tpm = db.Column(db.Numeric(precision=4, scale=1), nullable=False)


class Targets_summary(db.Model):
    isomir = db.Column(db.Text, db.ForeignKey('molecules.molecule', ondelete='CASCADE'), primary_key=True, index=True)
    cancer = db.Column(db.Text, db.ForeignKey('cancers.cancer', ondelete='CASCADE'), primary_key=True, index=True)
    
    num_targets_expressed = db.Column(db.SmallInteger)
    num_targets_anticorrelated = db.Column(db.SmallInteger)
    
    isomir_median_tpm = db.Column(db.Numeric(precision=4, scale=1), nullable=False)
    highly_expressed = db.Column(db.Boolean, nullable=False)
