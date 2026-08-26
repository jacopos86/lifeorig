import hashlib

#
#   read file
#

def _file_sha256(input_file):
    sha256 = hashlib.sha256()
    with open(input_file, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha256.update(block)
    return sha256.hexdigest()