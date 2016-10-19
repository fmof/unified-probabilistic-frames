import sys
load_minsky_attempts = 0
minsky_loaded = False
while load_minsky_attempts <= 1:
    try:
        import minsky
        minsky_loaded = True
        load_minsky_attempts += 1
    except Exception as e:
        if load_minsky_attempts == 0:
            import os
            script_dir = os.path.dirname(os.path.realpath(__file__))
            apath = os.path.realpath(os.path.join(script_dir, '../minsky/python'))
            sys.stderr.write("Could not find minsky module on path, trying again after adding\n\t%s\n" % (apath))
            sys.path.append(apath)
        load_minsky_attempts += 1
if not minsky_loaded:
    raise Exception("Could not find minsky module; path is: %s" % str(sys.path))
