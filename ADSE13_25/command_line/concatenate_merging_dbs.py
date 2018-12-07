from libtbx.phil import parse
import sys, os, shutil
from six.moves import StringIO

phil_scope = parse("""
  src_tag = None
    .type = strings
  src_dir = None
    .type = strings

  dest_tag = None
    .type = str
  dest_dir = "."
    .type = str
""")

def get_last_line(filename):
  # shamelessly lifted from stackoverflow
  with open(filename, "rb") as f:
      first = f.readline()        # Read the first line.
      f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
      while f.read(1) != b"\n":   # Until EOL is found...
	  f.seek(-2, os.SEEK_CUR) # ...jump back the read byte plus one more.
      last = f.readline()         # Read last line.
  return last

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception, e:
      print "Couldn't parse argument %s"%arg
      return

  scope = phil_scope.fetch(sources = user_phil)
  scope.show()
  params = scope.extract()

  dest_frames_db_filename = os.path.join(params.dest_dir, params.dest_tag + "_frame.db")
  dest_miller_db_filename = os.path.join(params.dest_dir, params.dest_tag + "_miller.db")
  dest_observation_db_filename = os.path.join(params.dest_dir, params.dest_tag + "_observation.db")

  for src_dir, src_tag in zip(params.src_dir, params.src_tag):
    source_frames_db_filename = os.path.join(src_dir, src_tag + "_frame.db")
    source_miller_db_filename = os.path.join(src_dir, src_tag + "_miller.db")
    source_observation_db_filename = os.path.join(src_dir, src_tag + "_observation.db")

    for filename in source_frames_db_filename, source_miller_db_filename, source_miller_db_filename:
      assert os.path.exists(filename), filename

    testing_filenames = [os.path.exists(f) for f in dest_frames_db_filename, dest_miller_db_filename, dest_observation_db_filename]
    assert testing_filenames.count(True) in [0,3]

    if testing_filenames.count(True) == 0:
      shutil.copy(source_miller_db_filename, dest_miller_db_filename)
      shutil.copy(source_frames_db_filename, dest_frames_db_filename)
      shutil.copy(source_observation_db_filename, dest_observation_db_filename)
      print "Copying %s to %s"%(src_tag, params.dest_tag)
      continue
    print "Appending %s to %s"%(src_tag, params.dest_tag)

    for line1, line2 in zip(open(source_miller_db_filename).readlines(), open(dest_miller_db_filename).readlines()):
      assert line1 == line2

    last_dest_frame = get_last_line(dest_frames_db_filename)
    n_dest_frames = int(last_dest_frame.split()[0])+1

    buf = StringIO()
    with open(source_frames_db_filename) as source_file:
      for line in source_file:
	data = line.strip().split()
	new_frame_id = str(int(data[0]) + n_dest_frames)
	buf.write(" ".join([new_frame_id] + data[1:]) + "\n")
    with open(dest_frames_db_filename, 'a') as dest_file:
      buf.seek (0)
      shutil.copyfileobj (buf, dest_file)

    buf = StringIO()
    with open(source_observation_db_filename) as source_file:
      for line in source_file:
	data = line.strip().split()
	new_frame_id = str(int(data[5]) + n_dest_frames)
	data[5] = new_frame_id
	buf.write(" ".join(data) + "\n")
    with open(dest_observation_db_filename, 'a') as dest_file:
      buf.seek (0)
      shutil.copyfileobj (buf, dest_file)

  print "Did it with style"

if __name__ == "__main__":
  run(sys.argv[1:])
