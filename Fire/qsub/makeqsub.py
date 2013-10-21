# Make commands that look like:
#   qsub -q daglab -l nodes=1:ppn=8 qsub_hivset_001.sh
#   001 -- 119

# script files look like:
#   matlab -nodisplay -r "cd /vision/u/liuyipei/joni/phylo/Fire/qsub/uri_flu/; matlabpool local 8; qsub_hivset_001;"

import subprocess
import os
import re

#batch = 'Uri0806fasta_90id'
batch = re.sub("^.*/", '', os.getcwd()).strip() # this returns the current directory
fasta_dir =           '/filesystem/u/liuyipei/vdj_data/%s/' % batch
Fire_qsub_batch_dir = '/filesystem/u/liuyipei/joni/phylo/Fire/qsub/%s' % batch
if re.match('/vision/', os.getcwd()):
  fasta_dir = re.sub('^/filesystem/', '/vision/', fasta_dir)
  Fire_qsub_batch_dir = re.sub('^/filesystem/', '/vision/', Fire_qsub_batch_dir)
if re.match('/visionnfs/', os.getcwd()):
  fasta_dir = re.sub('^/filesystem/', '/visionnfs/', fasta_dir)
  Fire_qsub_batch_dir = re.sub('^/filesystem/', '/visionnfs/', Fire_qsub_batch_dir)
if re.match('/scail/', os.getcwd()):
  fasta_dir = re.sub('^/filesystem/', '/scail/', fasta_dir)
  Fire_qsub_batch_dir = re.sub('^/filesystem/', '/scail/', Fire_qsub_batch_dir)

print "input fasta dir  :%s" % fasta_dir
print "qsub scripts dir: %s" % Fire_qsub_batch_dir

filenames = os.listdir(fasta_dir)
filenames = sorted(filter(lambda x: re.search('^.*\.fa(sta)*$', x), filenames)) 
num_fasta = len(filenames)

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

#num_fasta = 30  # for debugging this script
max_tiny = 10 # 5 sequences
max_small = 200 # 100 sequences
max_med = 2000 # 1000 sequences
filelengths = [file_len(fasta_dir + x) for x in filenames] # filenames are sorted alphabetically
parallel = [1 for x in filelengths]

sh_tiny_file = open('submit_t_%d.sh' % max_tiny, 'w+');
sh_small_file = open('submit_s_%d.sh' % max_small, 'w+');
sh_medium_file = open('submit_m_%d.sh'% max_med, 'w+');
sh_large_file = open('submit_larger_%d.sh' % max_med, 'w+');
for i in range(num_fasta):
  print(filelengths[i])

  curr_out = sh_tiny_file
  parallel[i] = 1
  if filelengths[i] > max_tiny:
    curr_out = sh_small_file
    parallel[i] = 1
  if filelengths[i] > max_small:
    curr_out = sh_medium_file
    parallel[i] = 4
  if filelengths[i] > max_med:
    curr_out = sh_large_file
    parallel[i] = 12 

  if parallel[i] > 1:
    curr_out.write("qsub -q daglab -l nodes=1:ppn=%d qsub_%s_%04.d.sh\n" % (parallel[i], batch, i))
  else:
    curr_out.write("qsub -q daglab qsub_%s_%04.d.sh\n" % (batch, i))
sh_tiny_file.close()
sh_small_file.close()
sh_medium_file.close()
sh_large_file.close()

master_mfile_name = "qsub_generic.m"

for i in range(num_fasta): 

  if filelengths[i]<= max_tiny:
    continue # too  many tiny files to be bothered with

  m_name = "qsub_%s_%04.d.m" % (batch, i)
  m_file = open(m_name, 'w+')
  sh_name = "qsub_%s_%04.d.sh" % (batch, i)
  sh_file = open(sh_name, 'w+')
  
  m_file.write("fasta_dir = '%s'\n" % fasta_dir)
  m_file.write("filename = '%s'\n" % filenames[i])
  m_file.write("fasta_file_line_count = %d\n" % filelengths[i])
  in_template_file = open(master_mfile_name, 'r')
  for line in in_template_file:
    m_file.write(line)

  if parallel[i] > 1:
    parallel_str = 'try matlabpool open local %d; catch \'unable to open local %d\', end;' % (parallel[i], parallel[i])
    end_parallel_str = 'try matlabpool close; catch \'unable to open close pool\', end;'
    sh_file.write('matlab_r2011b -nodisplay -nosplash -r "cd %s/; %s qsub_%s_%04.d;%s"' 
      % (Fire_qsub_batch_dir, parallel_str, batch, i, end_parallel_str))
  else:
    sh_file.write('matlab_r2011b -nodisplay -r "cd %s/; qsub_%s_%04.d;"' 
      % (Fire_qsub_batch_dir, batch, i))

  in_template_file.close()
  m_file.close()
  sh_file.close()
