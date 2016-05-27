#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

"""
###############################################################################
BamAlignmentViewer
A program for plot the reads covered target position for BAM sequence alignment
 files.
Author: JoeyHwong (joeyhwong@hotmail.com)
Licence: MIT
usage: plotreads [-h] [--version] --bam [bamfile] --region [region] [options]
positional arguments:
  -b, --bam     [FILE]  bam files containing mapped reads
  -r, --region  [STR]   the target region (0 base), example: chr1:123456-123457
optional arguments:
  -h, --help          show this help message and exit
  --version           show program's version number and exit
  --prefix    [STR]   string to add to the start of output graph.
  --rg        [STR]   The read group you want to draw
  -f          [NUM]   only show reads which mapping quality > NUM
  -d                  remove duplication reads, cigar level dups
  --ref       [FILE]  hg19 reference fasta file path
  -R                  show with reference, depend on -ref option specified

###############################################################################
"""

import sys
import os
import pysam
import re
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw
from optparse import OptionParser


class ReadPlot(object):
	def __init__(self, bam, outdir=os.getcwd(), readsgroup=None, qual_filter=1, rmdup=False, reference=None, **kwargs):
		self.outdir = outdir
		self.bam = self.prepare_index(bam, 'bam')
		self.reference = self.prepare_index(reference, 'fasta')
		self._bam = pysam.AlignmentFile(bam, "rb")
		self.Rmdup = rmdup
		self.qua_ctr = int(qual_filter)
		self.readsgroup = readsgroup
		self._refer = pysam.FastaFile(self.reference) if self.reference else None
		self.fontsize = int(kwargs["fontsize"]) if "fontsize" in kwargs else 10
		self.font = ImageFont.truetype(kwargs["font"], self.fontsize) if "font" in kwargs else ImageFont.load_default()
		self.leftpadding = int(kwargs["leftpadding"]) if "leftpadding" in kwargs else 3
		self.bgcolor = kwargs["bgcolor"] if "bgcolor" in kwargs else "#FFF"

	def __del__(self):
		self._bam.close()
		if self._refer:
			self._refer.close()

	def get_png(self, region, prefix):
		chrom, tmp_r = region.split(":")
		start, stop = map(int, tmp_r.split("-"))
		reads = self.get_reads(chrom, start, stop)
		f_out = "_".join(map(str, [prefix, chrom, start, stop])) + '.png'
		begin = max(start - 89, 1)
		line_weight, line_height = self.font.getsize(reads[0])[:2]
		img_height = line_height * (len(reads) + 4)
		gap = line_weight / len(reads[0])
		img_weight = line_weight + self.leftpadding * gap * 2
		img = Image.new("RGBA", (img_weight, img_height), self.bgcolor)
		draw = ImageDraw.Draw(img)
		y = line_height * 2
		x = (start - begin + self.leftpadding) * gap
		draw.rectangle(((gap, line_height), (img_weight - gap, img_height - line_height)), fill=None, outline="black")
		draw.line((x, y, x, y + line_height * len(reads)), fill="orange", width=(stop - start) * gap)
		x = self.leftpadding * gap
		for read in reads:
			if read[:2].isdigit():
				color = self.get_color(read[:2])
				read = read[2:]
			else:
				color = "#000"
			draw.text((x, y), read, color, font=self.font)
			y += line_height

		img.save(os.path.join(self.outdir, f_out))

	@staticmethod
	def get_color(number):
		code = hex(max(60 - int(number), 0))[2:]
		return "#" + code * 3

	def get_reads(self, chrom, start, stop):
		begin = max(start - 89, 1)
		end = stop + 89
		if stop == start:
			start -= 1
		if self.reference:
			end = min(end, self._refer.get_reference_length(chrom))
			refer = self._refer.fetch(chrom, begin, end)
		else:
			refer = "*" * (end - begin + 1)
		total_reads = list()
		indel = self.check_indel(chrom, start, stop)
		refer_bases = ""
		for i in range(begin, end):
			refer_bases += refer[i - begin]
			if i in indel and indel[i] > 0:
				refer_bases += "." * indel[i]
		total_reads.append(refer_bases)
		for aln_line in self._bam.fetch(chrom, start, stop):
			if self._filter_reads(aln_line) is True:
				continue
			cigar, pos_start = aln_line.cigartuples, aln_line.reference_start
			stand = "-" if aln_line.is_reverse is True else "+"
			sequece = aln_line.query_sequence
			base_qual = str(aln_line.mapping_quality).zfill(2)
			alignment_start = aln_line.query_alignment_start
			pos_start -= alignment_start
			reads = base_qual + " " * (pos_start - begin)
			pos = pos_start
			dele = 0
			for cigar_terms in cigar:
				if cigar_terms[0] == 0:
					base, number = 0, dele
					while base < cigar_terms[1]:
						if sequece[pos + base - pos_start] == refer_bases[pos - begin + number]:
							reads += str(stand)
							base += 1
						elif refer_bases[pos - begin + number] == ".":
							reads += "."
						else:
							reads += sequece[pos + base - pos_start]
							base += 1
						number += 1
					pos += cigar_terms[1]
				elif cigar_terms[0] == 1:
					for base in range(cigar_terms[1]):
						reads += sequece[pos - pos_start + base]
					pos += cigar_terms[1]
				elif cigar_terms[0] == 2:
					for base in range(cigar_terms[1]):
						reads += "^"
					dele += cigar_terms[1]
				else:
					for base in range(cigar_terms[1]):
						reads += sequece[pos - pos_start + base]
					pos += cigar_terms[1]
			total_reads.append(reads)
		return total_reads

	def _filter_reads(self, reads):
		if self.Rmdup and reads.is_duplicate:
			return True
		if reads.mapping_quality < self.qua_ctr or reads.is_qcfail:
			return True
		if self.readsgroup is not None and reads.get_tag('RG') != self.readsgroup:
			return True
		return False

	def check_indel(self, chrom, start, stop):
		indel = dict()
		for aln_line in self._bam.pileup(chrom, start, stop):
			indel[aln_line.reference_pos] = 0
			for pileup_read in aln_line.pileups:
				if self._filter_reads(pileup_read.alignment) is True:
					continue
				if pileup_read.indel > 0:
					indel[aln_line.reference_pos] = max(indel[aln_line.reference_pos], pileup_read.indel)
				else:
					continue
		return indel

	def prepare_index(self, in_file, preset='bam'):
		if not in_file or not os.path.isfile(in_file):
			return None
		if preset is 'bam':
			suffix = '.bai'
		elif preset is 'fasta':
			suffix = '.fai'
		else:
			suffix = '.tbi'
		if os.path.isfile(in_file + suffix) or os.path.isfile(re.sub("\.%s$" % preset, suffix, in_file)):
			return in_file
		index_dir, file_name = os.path.split(in_file)
		if not os.access(index_dir, os.W_OK):
			new_file = os.path.join(self.outdir, file_name)
			if os.path.exists():
				os.remove(new_file)
			os.symlink(in_file, new_file)
			in_file = new_file
		if preset is 'bam':
			pysam.index(in_file)
		elif preset is 'fasta':
			reffile = pysam.FastaFile(in_file)
			reffile.close()
		else:
			raise IndexError("Sorry, cann't load %s 's index ! " % in_file)
		return in_file


def main():
	usage = "Usage: %prog --bam [bamfile] --region [region] [options]"
	description = "A program for plot the reads covered target position for BAM sequence alignment files."
	author = "Author: joeyhwong@hotmail.com"
	parser = OptionParser(usage=usage, version="%prog 0.1", description=description,
	                      epilog=author)
	parser.add_option('-b', '--bam', help='The sorted bam file (required)',
	                  dest='bamfile')
	parser.add_option('-r', '--region', help='The target region (1 base), example: '
	                                         'chr1:123456-123457', dest='region')
	parser.add_option('--prefix', help='String to add to the start of output graph',
	                  default=os.path.join(os.getcwd(), 'Test'), dest='prefix')
	parser.add_option('--rg', help='The read group you want to draw', default=None,
	                  dest='readsgroup')
	parser.add_option('-f', help='Only show reads which mapping quality > NUM',
	                  default=1, type='int', dest='quality')
	parser.add_option('-d', help='Remove duplication reads, cigar level dups',
	                  action="store_true", dest="Rmdup", default=False)
	parser.add_option('--ref', help='The hg19 reference fasta file path', dest="refer", default=None)
	(options, args) = parser.parse_args()

	if not options.bamfile or not options.region:
		parser.print_help()
		return 'Expected arguments lost !!!'
	bamfile = os.path.abspath(options.bamfile)
	regions = str(options.region)
	assert os.path.isfile(bamfile), "%s not exit !" % bamfile
	prefix = os.path.abspath(options.prefix).strip('_')
	outdir, prefix = os.path.split(prefix)
	rg = options.readsgroup
	qua_ctl = int(options.quality)
	rmdup = options.Rmdup
	reference = os.path.abspath(options.refer) if options.refer else None
	test = ReadPlot(bamfile, outdir, rg, qua_ctl, rmdup, reference)
	test.get_png(region=regions, prefix=prefix)


if __name__ == '__main__':
	sys.exit(main())
