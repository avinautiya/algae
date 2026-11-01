"""ice_cliff_survival: standalone mathematical/executable foundations for stream-fed ice-cliff survival and reburial (milestone 0).
Everything here is a labelled structural baseline or a verified component; nothing is an environmental forecast. Unsupported runs raise UnsupportedRun."""
VERSION = "0.1.0-m0"
class UnsupportedRun(RuntimeError):
    """Raised when a requested run needs physics that is not implemented/verified (no empty success-shaped output)."""
