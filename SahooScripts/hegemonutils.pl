#!/usr/bin/perl -I /booleanfs/sahoo/scripts

if (scalar(@ARGV) <= 0) {
    print "perl convert.pl <cmd> <args>\n";
    exit(1);
}

use U;
use Hegemon;
use File::Basename;

my $cmd = shift(@ARGV);
if ($cmd eq "topcl") {
  &convertToPCL(@ARGV);
}
if ($cmd eq "toexpr") {
  &convertToExpr(@ARGV);
}
if ($cmd eq "toidx") {
  &convertToIdx(@ARGV);
}
if ($cmd eq "toih") {
  &convertToIh(@ARGV);
}
if ($cmd eq "tosurv") {
  &convertToSurv(@ARGV);
}
if ($cmd eq "Info") {
  &convertInfo(@ARGV);
}
if ($cmd eq "VInfo") {
  &convertVInfo(@ARGV);
}
if ($cmd eq "Surv") {
  &convertSurv(@ARGV);
}
if ($cmd eq "bv") {
  &convertBv(@ARGV);
}
if ($cmd eq "bvThr") {
  &convertBvThr(@ARGV);
}
if ($cmd eq "bs") {
  &analyzeBoolean(@ARGV);
}
if ($cmd eq "plot") {
  &analyzePlot(@ARGV);
}
if ($cmd eq "html") {
  &printHtml(@ARGV);
}
if ($cmd eq "print.mean") {
  &printAllTmpMean(@ARGV);
}
if ($cmd eq "tmp.mean") {
  &printTmpMean(@ARGV);
}
if ($cmd eq "tmp.files") {
  &printTmpFiles(@ARGV);
}
if ($cmd eq "gunzipAGCC") {
  &gunzipAGCC(@ARGV);
}
if ($cmd eq "orderArrays") {
  &orderArrays(@ARGV);
}
if ($cmd eq "corr") {
  &printCorrelation(@ARGV);
}
if ($cmd eq "corr1") {
  &printCorrelation1(@ARGV);
}
if ($cmd eq "gsea") {
  &printGsea(@ARGV);
}
if ($cmd eq "detect") {
  &printDetect(@ARGV);
}
if ($cmd eq "detectBv") {
  &printDetectBv(@ARGV);
}
if ($cmd eq "data") {
  &printData(@ARGV);
}
if ($cmd eq "seq2pcl") {
  &printSeqPCL(@ARGV);
}
if ($cmd eq "pcl2chip") {
  &printPCLChip(@ARGV);
}
if ($cmd eq "gsminfo") {
  &printGSMinfo(@ARGV);
}
if ($cmd eq "sort") {
  &printSort(@ARGV);
}
if ($cmd eq "shuffle") {
  &shuffleFile(@ARGV);
}
if ($cmd eq "heatmap") {
  &analyzeHeatmap(@ARGV);
}
if ($cmd eq "pair-info") {
  &analyzePairInfo(@ARGV);
}
if ($cmd eq "filter") {
  &analyzeFilter(@ARGV);
}
if ($cmd eq "merge") {
  &analyzeMerge(@ARGV);
}
if ($cmd eq "hist") {
  &analyzeHistogram(@ARGV);
}
if ($cmd eq "hidr") {
  &analyzeHighDR(@ARGV);
}
if ($cmd eq "hivdr") {
  &analyzeHighDynRange(@ARGV);
}
if ($cmd eq "hisdr") {
  &analyzeHighSDR(@ARGV);
}
if ($cmd eq "cluster") {
  &analyzeCluster(@ARGV);
}
if ($cmd eq "bs1") {
  &analyzeBoolean1(@ARGV);
}
if ($cmd eq "convert-xls") {
  &convertXLS(@ARGV);
}
if ($cmd eq "probe-info") {
  &printProbeInfo(@ARGV);
}
if ($cmd eq "count") {
    &analyzeCount(@ARGV);
}
if ($cmd eq "tpm") {
  &analyzeTPM(@ARGV);
}
if ($cmd eq "tpml") {
  &analyzeTPML(@ARGV);
}
if ($cmd eq "convert-geo") {
  &convertGEO(@ARGV);
}

sub convertToPCL {
  my $file = shift;
  &U::convertpcl($file, 2);
}

sub convertToExpr {
  my $file = shift;
  &U::convertexpr($file);
}

sub convertToIdx {
  my $file = shift;
  &U::convertidx($file);
}

sub convertToIh {
  my $file = shift;
  &U::convertih($file);
}

sub convertToSurv {
  my $file = shift;
  &U::convertsurv($file);
}

sub convertInfo {
  my ($dataset, $thr, @rest) = @_;
  my $pre = $dataset;
  $thr = "$pre\-thr.txt" if (!defined $thr);
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => $thr,
      survival => undef);
  $h->printInfo();
}

sub convertVInfo {
  my ($dataset, $thr, @rest) = @_;
  my $pre = $dataset;
  $thr = "$pre\-thr.txt" if (!defined $thr);
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => $thr,
      survival => undef);
  $h->printVInfo();
}

sub convertSurv {
  my ($dataset, $ct, $maxt, @rest) = @_;
  my $pre = $dataset;
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => "$pre\-survival.txt");
  $h->printAllSurvival($ct, $maxt);
}

sub convertBv {
  my ($dataset, @rest) = @_;
  my $pre = $dataset;
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => undef);
  $h->printBvFile();
}

sub getBvCode {
  my ($val, $thr, $athr) = @_;
  if (!defined $val || $val eq "") {
    return " ";
  }
  if ($val < $athr) {
    return " ";
  }
  elsif ($val < $thr->[2]) {
    return "0";
  }
  elsif ($val >= $thr->[3]) {
    return "2";
  }
  else {
    return "1";
  }
}

sub convertBvThr {
  my ($dataset, $athr, @rest) = @_;
  my $pre = $dataset;
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => undef);
  my $fh = $h->{'fh'};
  seek($fh, 0, 0);
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @header = split("\t", $head);
  my $start = $h->getStart();
  my $end = $h->getEnd();
  print join("\t", "ArrayID", "Name", "BitVector"), "\n";
  my $res = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $id = $list[0];
    my $name = $list[1];
    my $thr = $h->getThrData($id);
    my $bv = join("", map{&getBvCode($list[$_], $thr, $athr)} $start .. $end);
    print join("\t", $id, $name, $bv), "\n";
  }
}

sub analyzeBoolean {
  my ($pre, $n, $num, $index, $st, $p, $id3, @params) = @_;
  my $hithr = 50;
  my $sdthr = 0.5;
  my $drthr = 4;
  my $vdrthr;
  my ($thrx0, $thrx2, $thry0, $thry2);
  my $thrFile = "$pre\-thr.txt";
  my $infoFile = "$pre\-info.txt";
  my $vinfoFile = "$pre\-vinfo.txt";
  foreach (@params) {
    my ($k, $v) = split("=");
    if ($k eq "hithr") { $hithr = $v; }
    if ($k eq "sdthr") { $sdthr = $v; }
    if ($k eq "drthr") { $drthr = $v; }
    if ($k eq "vdrthr") { $vdrthr = $v; }
    if ($k eq "thrx0") { $thrx0 = $v; }
    if ($k eq "thrx2") { $thrx2 = $v; }
    if ($k eq "thry0") { $thry0 = $v; }
    if ($k eq "thry2") { $thry2 = $v; }
    if ($k eq "thrfile") { $thrFile = $v; }
    if ($k eq "infofile") { $infoFile = $v; }
    if ($k eq "vinfofile") { $vinfoFile = $v; }
  }
  my @rest = ($thrx0, $thrx2, $thry0, $thry2);
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => $thrFile,
      survival => undef);
  my $pGroups = undef;
  my $ihash;
  my $ihh = {};
  if (-e $infoFile) {
    $ihash = &U::getHash($infoFile);
    my @headers = @{$ihash->{'AffyID'}};
    for (my $i = 0; $i < scalar(@headers); $i++) {
      $ihh->{$headers[$i]} = $i;
    }
  }
  my $vihash;
  my $vihh = {};
  if (-e $vinfoFile) {
    $vihash = &U::getHash($vinfoFile);
    my @headers = @{$vihash->{'ProbeID'}};
    for (my $i = 0; $i < scalar(@headers); $i++) {
      $vihh->{$headers[$i]} = $i;
    }
  }
  
  my $l = $h->getIDs($n);
  for (my $i = 0; $i < scalar(@$l); $i++) {
    print STDERR "$n\[$i] = ", $l->[$i], "\n";
  }
  my $id = $l->[0];
  if ($num < 0) {
    $id = $h->getBestId($n);
  }
  else {
    $id = $l->[$num];
  }
  if (defined $id3 && $id3 ne '-' && !-e $id3) {
    $id3 = $h->getBestId($id3);
    my $bs = $h->getBooleanStats($id, $id3, $pGroups, @rest);
    my $n3 = $h->getName($id3);
    print join("\t", $id3, $n3, map { sprintf("%.2f", $_) } map { ($bs->[0]->[$_], $bs->[1]->[$_]) } 0 .. 3), "\n";
    print join("\t", $id3, $n3, map { sprintf("%.2f", $_) } map { ($bs->[2]->[$_], $bs->[3]->[$_]) } 0 .. 3), "\n";
    print "DynR = [id,name,thr,mean,mean-thr,perc,min,max,sd,thrNum,hi,int,lo]\n";
    my $res = $h->getDynamicRangeInfo($id3, $pGroups);
    print "DynR = [", join(",", $id3, $n3, (map { sprintf("%.2f", $res->[0]->[$_]) } 2 .. 8), map { $res->[0]->[$_] } 9 .. 12), "]\n";
    my $res = $h->getDynamicRangeInfo($id, $pGroups);
    print "DynR = [", join(",", $id, $n, (map { sprintf("%.2f", $res->[0]->[$_]) } 2 .. 8), map { $res->[0]->[$_] } 9 .. 12), "]\n";
    if (defined $vihash) {
      my $vdr3 = $vihash->{$id3}->[$vihh->{'Perc 0.95'}] - $vihash->{$id3}->[$vihh->{'Perc -0.95'}];
      my $sdr1 = $vihash->{$id3}->[$vihh->{'Perc 0.95'}] - $vihash->{$id3}->[$vihh->{'thr'}];
      my $sdr2 = $vihash->{$id3}->[$vihh->{'thr'}] - $vihash->{$id3}->[$vihh->{'Perc -0.95'}];
      my $s3 = abs($sdr1/$sdr2);
      print "vDynR = [$id3 = $vdr3, $s3, $sdr1, $sdr2]\n";

      my $vdr1 = $vihash->{$id}->[$vihh->{'Perc 0.95'}] - $vihash->{$id}->[$vihh->{'Perc -0.95'}];
      my $sdr1 = $vihash->{$id}->[$vihh->{'Perc 0.95'}] - $vihash->{$id}->[$vihh->{'thr'}];
      my $sdr2 = $vihash->{$id}->[$vihh->{'thr'}] - $vihash->{$id}->[$vihh->{'Perc -0.95'}];
      my $s1 = abs($sdr1/$sdr2);
      print "vDynR = [$id = $vdr1, $s1, $sdr1, $sdr2]\n";
    }
    return;
  }
  my @ids = keys(%{$h->{"idHash"}});
  if (-e $id3) {
    my $hh = &U::getHash($id3);
    @ids = keys(%{$hh});
  }
  print join("\t", "AID", "Name", map { ("ST[$_]", "p[$_]") } 0 .. 3), "\n";
  print join("\t", "$id", "$n\[$num]", map { ("", "") } 0 .. 3), "\n";
  my $idx = 0;
  for (my $j = 0; $j < scalar(@ids) ; $j++) {
    my $id2 = $ids[$j];
    next if (!defined $h->{"idHash"}->{$id2});
    if (defined $ihash) {
      next if ($ihash->{$id2}->[$ihh->{'hi'}] < $hithr);
      next if ($ihash->{$id2}->[$ihh->{'sd'}] < $sdthr);
      my $dr = $ihash->{$id2}->[$ihh->{'max'}] - $ihash->{$id2}->[$ihh->{'min'}];
      next if ($dr < $drthr);
    }
    if (defined $vihash) {
      my $vdr = $vihash->{$id2}->[$vihh->{'Perc 0.95'}] - $vihash->{$id2}->[$vihh->{'Perc -0.95'}];
      next if (defined $vdrthr && $vdr < $vdrthr);
    }
    my $n2 = $h->getName($id2);
    my $bs = $h->getBooleanStats($id, $id2, $pGroups, @rest);
    my $str = join("\t", $id2, $n2, map { sprintf("%.2f", $_) } 
        map { ($bs->[2]->[$_], $bs->[3]->[$_]) } 0 .. 3);
    
    if (defined $index && defined $st && defined $p) {
      if ($index <= 3 && $bs->[2]->[$index] > $st && $bs->[3]->[$index] <= $p) {
        print $str, "\n";
      }
      elsif ($index == 4 && $bs->[2]->[1] > $st && $bs->[3]->[1] <= $p
          && $bs->[2]->[2] > $st && $bs->[3]->[2] <= $p) {
        print $str, "\n";
      }
      elsif ($index == 5 && $bs->[2]->[0] > $st && $bs->[3]->[0] <= $p
          && $bs->[2]->[3] > $st && $bs->[3]->[3] <= $p) {
        print $str, "\n";
      }
    }
    else {
      print $str, "\n";
    }
    if (($idx % 1000) == 0) {
        print STDERR "$idx\n";
    }
    $idx++;
  }
}

sub analyzePlot {
  my ($pre, $n1, $n2, $num1, $num2, $outfile, $thr1, $thr2, @params) = @_;
  $num1 = -1 if (!defined $num1);
  $num2 = -1 if (!defined $num2);
  $outfile = "plot.png" if (!defined $outfile);
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => undef);
  my $l = $h->getIDs($n1);
  for (my $i = 0; $i < scalar(@$l); $i++) {
    print STDERR "$n1\[$i] = ", $l->[$i], "\n";
  }
  my $id1 = $l->[0];
  if ($num1 < 0) {
    $id1 = $h->getBestId($n1);
    print STDERR "using $n1\[$num1] = ", $id1, "\n";
  }
  else {
    $id1 = $l->[$num1];
    print STDERR "using $n1\[$num1] = ", $id1, "\n";
  }
  my $l = $h->getIDs($n2);
  for (my $i = 0; $i < scalar(@$l); $i++) {
    print STDERR "$n2\[$i] = ", $l->[$i], "\n";
  }
  my $id2 = $l->[0];
  if ($num2 < 0) {
    $id2 = $h->getBestId($n2);
    print STDERR "using $n2\[$num2] = ", $id2, "\n";
  }
  else {
    $id2 = $l->[$num2];
    print STDERR "using $n2\[$num2] = ", $id2, "\n";
  }
  $thr1 = [$thr1, 3, $thr1-0.5, $thr1+0.5] if (defined $thr1);
  $thr2 = [$thr2, 3, $thr2-0.5, $thr2+0.5] if (defined $thr2);
  $h->plotBooleanPair($outfile, $id1, $id2, undef, $thr1, $thr2);
}

sub printHtml {
  my ($cmd, $pre, $id, $rest) = @_;
  my $idlist = [];
  if (defined $id) {
    push @$idlist, $id;
  }
  while (<STDIN>) {
    s/[\r\n]//g;
    next if (/^ID/);
    if (/Using.*= (.*)/) {
      unshift @$idlist, $1;
      next;
    }
    my @list = split("\t");
    push @$idlist, $list[0];
  }
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => undef);
  $h->subplot($cmd, $idlist);
}

sub printAllTmpMean {
  my @files = @_;
  my @fhs;
  foreach my $i (0 .. $#files) {
    my $file = $files[$i];
    open(my $fh, "<$file") || die "Can't open $file\n";
    $fhs[$i] = $fh;
  }
  print "NumChips";
  my $total = 0;
  foreach my $i (0 .. $#files) {
    my $fh = $fhs[$i];
    my $buf;
    read($fh, $buf, 4);
    my $nc = unpack("i", $buf);
    $total += $nc;
    print "\t$nc";
  }
  print "\t$total\n";
  print "NumProbes";
  my $total = 0;
  foreach my $i (0 .. $#files) {
    my $fh = $fhs[$i];
    my $buf;
    read($fh, $buf, 4);
    my $np = unpack("i", $buf);
    $total += $np;
    print "\t$np";
  }
  print "\t$total\n";
  my $index = 0;
  while (1) {
    my $done = 1;
    my $total = 0;
    print "Mean[$index]";
    foreach my $i (0 .. $#files) {
      my $fh = $fhs[$i];
      my $buf;
      my $n = read($fh, $buf, 8);
      if ($n == 8) {
        my $val = unpack("d", $buf);
        $total += $val;
        print "\t$val";
        $done = 0;
      }
      else {
        print "\t";
      }
    }
    print "\t$total\n";
    $index++;
    last if ($done == 1);
  }
  foreach my $i (0 .. $#files) {
    my $fh = $fhs[$i];
    close($fh);
  }
}

sub printTmpMean {
  my @files = @_;
  my $numchips = 0;
  my $numprobes;
  my @means;
  foreach my $file (@files) {
    open(my $fh, "<$file") || die "Can't open $file\n";
    my $buf;
    read($fh, $buf, 4);
    $numchips += unpack("i", $buf);
    read($fh, $buf, 4);
    my $num = unpack("i", $buf);
    if (!defined $numprobes) {
      $numprobes = $num;
      @means = map { 0 } 0 .. ($numprobes - 1);
    }
    elsif ($numprobes != $num) {
      print STDERR "$file probes $num != $numprobes\n";
      exit(1);
    }
    for (my $i = 0; $i < $numprobes; $i++) {
      read($fh, $buf, 8);
      my $val = unpack("d", $buf);
      $means[$i] += $val;
    }
    close($fh);
  }
  print STDERR "$numchips $numprobes\n";
  print pack("i", $numchips), pack("i", $numprobes);
  for (my $i = 0; $i < $numprobes; $i++) {
    print pack("d", $means[$i]);
  }
}

sub gunzipAGCC {
  my ($file, $rest) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  while (<$fh>) {
    chomp;
    my @list = split("\t");
    if ($list[0] =~ /AGCC/ && $list[1] =~ /.gz$/) {
      system("gunzip $list[1] >/dev/null");
      $list[1] =~ s/.gz$//g;
    }
    print join("\t", @list), "\n";
  }
  close($fh);
}

sub orderArrays {
  my ($pclfile, $listfile, $rest) = @_;
  open(FL, "<$listfile") || die "Can't open $listfile\n";
  my @list;
  while (<FL>) {
    s/[\r\n]//g;
    my ($id, $rest) = split("\t");
    push @list, $id;
  }
  close(FL);

  open(FL, "<$pclfile") || die "Can't open $pclfile\n";
  my $head = <FL>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  my %hash;
  for (my $i = 0; $i < scalar(@headers); $i++) {
    $hash{$headers[$i]} = $i;
  }
  my $index = 0;
  foreach (@list) {
    if (defined $hash{$_}) {
      if ($index == 0) {
        print $headers[$hash{$_}];
        $index++;
      }
      else {
        print "\t", $headers[$hash{$_}];
      }
    }
    else {
      print STDERR "$_ not present\n";
    }
  }
  print "\n";
  while (<FL>) {
    s/[\r\n]//g;
    my @array = split("\t");
    my $first = 0;
    foreach (@list) {
      if (defined $hash{$_}) {
        if ($first == 0) {
          print $array[$hash{$_}];
          $first++;
        }
        else {
          print "\t", $array[$hash{$_}];
        }
      }
    }
    print "\n";
    $index++;
    if (($index % 1000) == 0) {
      print STDERR "$index\n";
    }
  }

  close(FL);
}

sub printCorrelation {
  my ($pre, $n, $num, @params) = @_;
  my $hithr = 50;
  my $sdthr = 0.5;
  my $drthr = 4;
  my $corrthr = 0.4;
  foreach (@params) {
    my ($k, $v) = split("=");
    if ($k eq "hithr") { $hithr = $v; }
    if ($k eq "sdthr") { $sdthr = $v; }
    if ($k eq "drthr") { $drthr = $v; }
    if ($k eq "corrthr") { $corrthr = $v; }
  }
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => undef);
  my $l = $h->getIDs($n);
  for (my $i = 0; $i < scalar(@$l); $i++) {
    print STDERR "$n\[$i] = ", $l->[$i], "\n";
  }
  my $id = $l->[0];
  if ($num < 0) {
    $id = $h->getBestId($n);
    print STDERR "using $n\[$num] = ", $id, "\n";
  }
  else {
    $id = $l->[$num];
    print STDERR "using $n\[$num] = ", $id, "\n";
  }
  my $group1 = $h->getArraysAll();
  my $cres = $h->correlation($id, $group1);
  my $chash = {};
  foreach my $i (@$cres) {
    $chash->{$i->[1]} = $i->[0];
  }
  if (-e "$pre\-info.txt") {
    my $ihash = &U::getHash("$pre\-info.txt");
    my @headers = @{$ihash->{'AffyID'}};
    my $ihh = {};
    for (my $i = 0; $i < scalar(@headers); $i++) {
      $ihh->{$headers[$i]} = $i;
    }
    foreach my $pi (keys %{$chash}) {
      my @n = split(" /// ", $h->getName($pi));
      my $n = $n[0];
      next if ($ihash->{$pi}->[$ihh->{'hi'}] < $hithr);
      next if ($ihash->{$pi}->[$ihh->{'sd'}] < $sdthr);
      my $dr = $ihash->{$pi}->[$ihh->{'max'}] - $ihash->{$pi}->[$ihh->{'min'}];
      next if ($dr < $drthr);
      next if ($chash->{$pi} < $corrthr && $chash->{$pi} >= -$corrthr);
      my $str = join("\t", $chash->{$pi}, $pi, $n, $dr);
      print $str, "\n";
    }
  }
  else {
    foreach my $i (@$cres) {
      print join("\t", @$i), "\n";
    }
  }
}

sub printCorrelation1 {
  my ($pre, $n, $num, @params) = @_;
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => undef);
  my $l = $h->getIDs($n);
  for (my $i = 0; $i < scalar(@$l); $i++) {
    print STDERR "$n\[$i] = ", $l->[$i], "\n";
  }
  my $id = $l->[0];
  if ($num < 0) {
    $id = $h->getBestId($n);
    print STDERR "using $n\[$num] = ", $id, "\n";
  }
  else {
    $id = $l->[$num];
    print STDERR "using $n\[$num] = ", $id, "\n";
  }
  open(my $fh, "<$pre\-idx.txt");
  my $count = 0;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    last if ($list[0] eq $id);
    $count++;
  }
  close($fh);
  my $cmd = "java -cp /booleanfs/sahoo/softwares/BooleanNet/dist/lib/tools.jar -Xms64m -Xmx4G tools.CustomAnalysis corrOne corr1234.txt $pre\-expr.txt $count > /dev/null";
  system("$cmd");
  open(my $fh, "<corr1234.txt");
  my $cres;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    push @$cres, [@list];
  }
  close($fh);
  foreach my $i (sort { $b->[0] <=> $a->[0] } @$cres) {
    print join("\t", @$i), "\n";
  }
}

sub printTmpFiles {
  my $listfile = shift;
  my $dir = &dirname($listfile);
  open(my $fh1, "<$listfile") || die "Can't open $listfile\n";
  while (<$fh1>) {
    s/[\r\n]//g;
    my $fdir = &dirname($_);
    my $file = &basename($_);
    if ($fdir =~ /^\./) {
      $fdir =~ s/^\./$dir/g;
    }
    open(my $fh, "<$fdir/$file\.mean") || die "Can't open $fdir/$file\.mean\n";
    my $buf;
    read($fh, $buf, 4);
    my $numchips = unpack("i", $buf);
    read($fh, $buf, 4);
    my $numprobes = unpack("i", $buf);
    close($fh);
    open(my $fh, "<$fdir/$file") || die "Can't open $fdir/$file\n";
    for (my $i = 0; $i < $numchips; $i++) {
      my $fs = 128;
      my $ns = $fs/4;
      my $ptr = (2 + $ns + $i * ($ns + $numprobes)) * 4;
      seek($fh, $ptr, 0);
      read($fh, $buf, $fs);
      my $val = unpack("Z*", $buf);
      print $val, "\n";
    }
  }
  close($fh1);
}

sub printGsea {
  my ($file, $prefix, $start, $end, $skip, $rest) = @_;
  $prefix = "expr" if (!defined $prefix);
  my $headers = &U::getHeaders($file);
  $start = 2 if (!defined $start);
  $end = scalar(@$headers) - 1 if (!defined $end);
  $skip = 0 if (!defined $skip);
  my $num = $end - $start + 1;

  my $cls = "$prefix\.cls";
  open(my $ocls, ">$cls") || die "Can't write $cls\n";
  print $ocls "$num 2 1\n";
  print $ocls "# 1 2\n";
  print $ocls join(" ", map { 1 } 0 .. ($num/2 - 1)), " ", 
        join(" ", map { 2 } 0 ..  ($num - $num/2 - 1)), "\n";
  close($ocls);

  my $chp = "$prefix\.chip";
  open(my $ochp, ">$chp") || die "Can't write $chp\n";
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  print $ochp "Probe Set ID\tGene Symbol\tGene Title\n";
  for (my $i = 0; $i < $skip; $i++) {
    my $h = <$fh>;
  }
  my $numprobes = 0;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $id = $list[0];
    my ($name, $desc) = split(/: /, $list[1]);
    my @n = split(" /// ", $name);
    print $ochp "$id\t$n[0]\t$desc\n";
    $numprobes++;
  }
  close($fh);
  close($ochp);

  my $rnk = "$prefix\.rnk";
  open(my $ornk, ">$rnk") || die "Can't write $rnk\n";
  my $gct = "$prefix\.gct";
  open(my $ogct, ">$gct") || die "Can't write $gct\n";
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  print $ogct "#1.2\n$numprobes $num\n";
  print $ogct "ProbeID\tName\t", 
        join("\t", map { $headers->[$_] } $start ..  $end), "\n";
  for (my $i = 0; $i < $skip; $i++) {
    my $h = <$fh>;
  }
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $id = $list[0];
    my $name = $list[1];
    print $ogct "$id\t$name\t",
        join("\t", map { $list[$_] } $start ..  $end), "\n";;
    my $m1 = &U::mean([map { $list[$_] } $start ..  ($start + $num/2 - 1)]);
    my $s = $start + $num/2;
    my $m2 = &U::mean([map { $list[$_] } $s ..  $end]);
    print $ornk "$id\t", ($m2 - $m1), "\n";
  }
  close($fh);
  close($ogct);
  close($ornk);
}

sub printDetect {
  my ($file, @files) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  my @fhs;
  foreach my $f (@files) {
    open(my $fh1, "<$f") || die "Can't open $f\n";
    push @fhs, $fh1;
  }
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  my @headers_target;
  foreach my $fh1 (@fhs) {
    my $head = <$fh1>;
    $head =~ s/[\r\n]//g;
    my @h = split("\t", $head);
    push @headers_target, @h;
  }

  my $hash = {};
  $hash->{"sum1"} = [map { 0 } @headers];
  $hash->{"sum2"} = [map { 0 } @headers_target];
  $hash->{"sum1sq"} = [map { 0 } @headers];
  $hash->{"sum2sq"} = [map { 0 } @headers_target];
  $hash->{"cov"} = [map { [map { 0 } @headers_target] } @headers];
  $hash->{"cor"} = [map { [map { 0 } @headers_target] } @headers];
  $hash->{"num"} = 0;
  while (<$fh>) {
    chomp;
    my @list1 = split("\t");
    my @list2;
    foreach my $fh1 (@fhs) {
      my $in = <$fh1>;
      chomp $in;
      my @h = split("\t", $in);
      push @list2, @h;
    }
    for (my $i = 0; $i < scalar(@headers); $i++) {
      $hash->{"sum1"}->[$i] += $list1[$i];
      $hash->{"sum1sq"}->[$i] += $list1[$i] * $list1[$i];
    }
    for (my $i = 0; $i < scalar(@headers_target); $i++) {
      $hash->{"sum2"}->[$i] += $list2[$i];
      $hash->{"sum2sq"}->[$i] += $list2[$i] * $list2[$i];
    }
    for (my $i = 0; $i < scalar(@headers); $i++) {
      for (my $j = 0; $j < scalar(@headers_target); $j++) {
        $hash->{"cov"}->[$i]->[$j] += $list1[$i] * $list2[$j];
      }
    }
    if ( ($hash->{"num"} % 100) == 0 ) {
      print STDERR $hash->{"num"}, "\n";
    }
    $hash->{"num"}++;
  }
  for (my $i = 0; $i < scalar(@headers); $i++) {
    my $max = 0;
    my $maxi = 0;
    for (my $j = 0; $j < scalar(@headers_target); $j++) {
      my $cov = $hash->{"cov"}->[$i]->[$j] * $hash->{"num"};
      $cov -= $hash->{"sum1"}->[$i] * $hash->{"sum2"}->[$j];
      my $n1 = ($hash->{"num"} * $hash->{"sum1sq"}->[$i]) - 
        ($hash->{"sum1"}->[$i]**2);
      my $n2 = ($hash->{"num"} * $hash->{"sum2sq"}->[$j]) - 
        ($hash->{"sum2"}->[$j]**2);
      my $n = $n1 * $n2;
      my $c = 0;
      if ($n > 0) { $c = $cov/sqrt($n); }
      $hash->{"cor"}->[$i]->[$j] = $c;
      if ($max < $c) { $max = $c; $maxi = $j; }
    }
    print "$headers[$i] -> $max $headers_target[$maxi]\n";
    my @sorted = sort { 
      $hash->{"cor"}->[$i]->[$b] <=> $hash->{"cor"}->[$i]->[$a]
    } 0 .. $#headers_target;
    foreach (0 .. 1000) {
      next if ($_ > $#headers_target);
      print join("\t", "", $hash->{"cor"}->[$i]->[$sorted[$_]], 
          $headers_target[$sorted[$_]]), "\n";
    }
  }
  close($fh);
  foreach my $fh1 (@fhs) {
    close($fh1);
  }
}

sub printDetectBv {
  my ($pclfile, $bvfile, @files) = @_;
  open(my $fh, "<$bvfile") || die "Can't open $bvfile\n";
  my @fhs;
  for (my $i = 1; $i < scalar(@files); $i+=2) {
    my $f = $files[$i];
    open(my $fh1, "<$f") || die "Can't open $f\n";
    push @fhs, $fh1;
  }
  my $h = &U::getHeaders($pclfile);
  my @headers = map { $h->[$_] } 3 .. (scalar(@$h) - 1);
  my @headers_target;
  for (my $i = 0; $i < scalar(@files); $i+=2) {
    my $h = &U::getHeaders($files[$i]);
    my @h = map { $h->[$_] } 3 .. (scalar(@$h) - 1);
    push @headers_target, @h;
  }

  my $hash = {};
  $hash->{"aub"} = [map { [map { 0 } @headers_target] } @headers];
  $hash->{"anb"} = [map { [map { 0 } @headers_target] } @headers];
  $hash->{"cor"} = [map { [map { 0 } @headers_target] } @headers];
  $hash->{"num"} = 0;
  while (<$fh>) {
    chomp;
    my @list1 = split("\t");
    my @list1 = split("", $list1[2]);
    my @list2;
    foreach my $fh1 (@fhs) {
      my $in = <$fh1>;
      chomp $in;
      my @h = split("\t", $in);
      my @h = split("", $h[2]);
      push @list2, @h;
    }
    for (my $i = 0; $i < scalar(@headers); $i++) {
      for (my $j = 0; $j < scalar(@headers_target); $j++) {
        if ($list1[$i] eq 2 && $list2[$j] eq 2) {
          $hash->{"anb"}->[$i]->[$j] ++;
        }
        if ($list1[$i] eq 2 || $list2[$j] eq 2) {
          $hash->{"aub"}->[$i]->[$j] ++;
        }
      }
    }
    if ( ($hash->{"num"} % 100) == 0 ) {
      print STDERR $hash->{"num"}, "\n";
    }
    $hash->{"num"}++;
  }
  for (my $i = 0; $i < scalar(@headers); $i++) {
    my $max = 0;
    my $maxi = 0;
    for (my $j = 0; $j < scalar(@headers_target); $j++) {
      my $n = $hash->{"anb"}->[$i]->[$j];
      my $d = $hash->{"aub"}->[$i]->[$j];
      my $c = 0;
      if ($d > 0) { $c = $n/$d; }
      $hash->{"cor"}->[$i]->[$j] = $c;
      if ($max < $c) { $max = $c; $maxi = $j; }
    }
    print "$headers[$i] -> $max $headers_target[$maxi]\n";
    my @sorted = sort { 
      $hash->{"cor"}->[$i]->[$b] <=> $hash->{"cor"}->[$i]->[$a]
    } 0 .. $#headers_target;
    foreach (0 .. 9) {
      next if ($_ > $#headers_target);
      print join("\t", "", $hash->{"cor"}->[$i]->[$sorted[$_]], 
          $headers_target[$sorted[$_]]), "\n";
    }
  }
  close($fh);
  foreach my $fh1 (@fhs) {
    close($fh1);
  }
}

sub printData {
  my ($dataset, @ids) = @_;
  my $pre = $dataset;
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => undef);
  my $start = $h->getStart();
  my $end = $h->getEnd();
  print join("\t", "len", map { $end - $start + 1 } @ids), "\n";
  my $thrData = [];
  foreach (@ids) {
    my $thr = $h->getThrData($_);
    push @$thrData, $thr;
  }
  my $names = [];
  foreach (@ids) {
    my $n = $h->getName($_);
    push @$names, "$_: $n";
  }
  print join("\t", "thr1", map { $_->[0] } @$thrData), "\n";
  print join("\t", "thr0", map { $_->[2] } @$thrData), "\n";
  print join("\t", "thr2", map { $_->[3] } @$thrData), "\n";
  print join("\t", "names", map { $_ } @$names), "\n";
  my $exprData = [];
  foreach (@ids) {
    my $e = $h->getExprData($_);
    push @$exprData, $e;
  }
  foreach my $i ($start .. $end) {
    print join("\t", $i - $start + 1, map { $_->[$i] } @$exprData), "\n";
  }
}

sub printSeqPCL {
  my $file = shift;
  open(my $fh, "<$file");
  my $header = <$fh>;
  $header =~ s/[\r\n]//g;
  my @hdrs = split("\t", $header);
  my @orders;
  foreach my $i (0 .. $#hdrs) {
    if ($hdrs[$i] =~ /_FPKM$/) {
        push @orders, $i;
    }
  }
  print join("\t", "Probe Set ID\tGene Symbol: Gene Title\tGWEIGHT",
        map { $hdrs[$_] } @orders), "\n";
    print join("\t", "EWEIGHT\t\t", map { 1 } @orders), "\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    print join("\t", $list[0], $list[3].": ".$list[4], 1, 
        map { $list[$_] } @orders), "\n";
  }
  close($fh);
}

sub printPCLChip {
  my $file = shift;
  open(my $fh, "<$file");
  my $header = <$fh>;
  print "Probe Set ID\tGene Symbol\tGene Title\n";
  my $ew = <$fh>;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my ($name, $desc) = split(/: /, $list[1], 2);
    print join("\t", $list[0], $name, $desc), "\n";
  }
  close($fh);
}

sub printGSMinfo {
  my ($file, $expr, $pre, @rest) = @_;
  my $hash = {};
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $id;
  my $headers = ["title", "source", "sampleType", "SampleID", "Tissue",
    "Condition", "Expt"];
  my $phash = {};
  foreach my $k (@$headers) {
    $phash->{$k} = 1;
  }
  while (<$fh>) {
    s/[\r\n]//g;
    if (/^\^SAMPLE = (.*)$/) {
      $id = $1;
    }
    if (/^!Sample_type = (.*)$/) {
      $hash->{$id}->{'sampleType'} = $1;
    }
    if (/^!Sample_title = (.*)$/) {
      $hash->{$id}->{'title'} = $1;
    }
    if (/^!Sample_series_id = (.*)$/) {
      $hash->{$id}->{'source'} = $1;
    }
    if (/^!Sample_characteristics_ch1 = (.*)$/) {
      my @l = split(";", $1);
      foreach (@l) {
        if (/:/) {
          my ($a, $b) = split(/: /, $_, 2);
          next if (length($a) > 50);
          next if ($a =~ /:/);
          my $ctype = lc($a);
          $hash->{$id}->{$ctype} = $b;
          if (!defined $phash->{$ctype}) {
            push @$headers, $ctype;
          }
          $phash->{$ctype} = 1;
        }
      }
    }
  }
  close($fh);
  my $hdr = &U::getHeaders($expr);
  my @keys = @$headers;

  my $hh = {};
  foreach my $i (0 .. $#keys) {
    $hh->{$keys[$i]} = $i;
  }

  if (scalar(@rest) > 0 && defined $rest[0] && -f $rest[0]) {
    my $hash2 = &U::getHash($rest[0]);
    foreach my $id (@$hdr) {
      if ($id =~ /(GSM[0-9]+)/i) {
        my $arr = uc($1);
        if (defined $hash2->{$arr}) {
          foreach my $k (0 .. $#{$hash2->{"GSMID"}}) {
            next if (!defined $phash->{$hash2->{"GSMID"}->[$k]});
            $hash->{$arr}->{$hash2->{"GSMID"}->[$k]} = $hash2->{$arr}->[$k];
          }
        }
      }
    }
  }

  open(my $fh1, ">$pre\-survival.txt");
  open(my $fh2, ">$pre\-ih.txt");
  print $fh1 join("\t", "ArrayID", "time", "status", map { "c $_" } @keys), "\n";
  print $fh2 "ArrayID\tArrayHeader\tClinicalhHeader\n";
  foreach my $id (@$hdr) {
    if ($id =~ /(GSM[0-9]+)/i) {
        print $fh1 join("\t", $id, "", "", map { $hash->{uc($1)}->{$_} } @keys), "\n";
        print $fh2 join("\t", $id, $id, $hash->{uc($1)}->{'title'}), "\n";
    }
  }
  close($fh1);
  close($fh2);
}

sub printSort {
  my ($file, $index, $index2, $rest) = @_;
  my $hash = &U::getHash($file, $index2);
  my @keys = sort { $hash->{$b}->[$index] <=> $hash->{$a}->[$index] } keys(%{$hash});
  foreach my $k (@keys) {
    print join("\t", @{$hash->{$k}}), "\n";
  }
}

sub shuffleFile {
  my ($pclfile, @rest) = @_;
  my $ifh;
  open($ifh, "<$pclfile") || die "Can't open $pclfile\n";
  my $index = 0;
  while ( $in = <$ifh> ) {
    if ( ($index % 100) == 0) {
        print STDERR "[$index]\n";
    }
    $in =~ s/[\r\n]//g;
    my @list = split("\t", $in, -1);
    my @perm = (0 .. 2, &U::getRandPerm(3, $#list));
    @list = map { $list[$_] } @perm;
    print join("\t", @list)."\n";
    $index++;
  }
  close($ofh);
}

sub analyzeHeatmap {
  my ($file, @rest) = @_;
  my $expr = [];
  open (my $fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @hdr = split("\t", $head);
  my $start = 1;
  my $end = $#hdr;
  if ($hdr[1] =~ /Name/i) {
    $start = 2;
  }
  if ($hdr[2] eq "GWEIGHT") {
    $start = 3;
    my $ew = <$fh>;
  }
  my $columns = [map { $hdr[$_] } $start .. $end];
  my $rows = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    if ($start eq 1) {
      push @$rows, $list[0];
    }
    else {
      my ($name, $desc) = split(/: /, $list[1], 2);
      push @$rows, $name;
    }
    push @$expr, [map { $list[$_] } $start .. $end];
  }
  close($fh);
  my %params;
  foreach my $i (@rest) {
    my ($k, $v) = split("=", $i, 2);
    $params{$k} = $v;
  }
  my $outfile = "heatmap.png";
  if (defined $params{"ofile"}) {
    $outfile = $params{"ofile"};
    delete $params{"ofile"};
  }
  my $ph = undef;
  if (defined $params{"ph"}) {
    $ph = [undef, \*STDOUT];
    delete $params{"ph"};
  }
  &Hegemon::heatmap1($outfile, $expr, $rows, $columns, $ph, %params);
}

sub analyzePairInfo {
  my ($pre, $id1, $id2, $rest) = @_;
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => undef,
      survival => undef);
  my $bid1 = $h->getBestId($id1);
  my $bid2 = $h->getBestId($id2);
  my $x = $h->getExprData($bid1);
  my $y = $h->getExprData($bid2);
  print "$id1: $bid1\n";
  print "$id2: $bid2\n";
  $x->[0] = ""; $x->[1] = "";
  $y->[0] = ""; $y->[1] = "";
  my $hash = &U::getXYstats($x, $y);
  foreach my $k (keys %{$hash}) {
    print join("\t", $k, $hash->{$k}), "\n";
  }
}

sub analyzeFilter {
  my ($pre, @params) = @_;
  my $hithr;
  my $sdthr;
  my $drthr;
  my $vdrthr;
  my $sdrthr;
  my $hidrthr;
  my $thrFile = "$pre\-thr.txt";
  my $infoFile = "$pre\-info.txt";
  my $vinfoFile = "$pre\-vinfo.txt";
  foreach (@params) {
    my ($k, $v) = split("=");
    if ($k eq "hithr") { $hithr = $v; }
    if ($k eq "sdthr") { $sdthr = $v; }
    if ($k eq "drthr") { $drthr = $v; }
    if ($k eq "vdrthr") { $vdrthr = $v; }
    if ($k eq "sdrthr") { $sdrthr = $v; }
    if ($k eq "hidrthr") { $hidrthr = $v; }
    if ($k eq "thrfile") { $thrFile = $v; }
    if ($k eq "infofile") { $infoFile = $v; }
    if ($k eq "vinfofile") { $vinfoFile = $v; }
  }
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => $thrFile,
      survival => undef);
  my $ihash;
  my $ihh = {};
  if (-e $infoFile) {
    $ihash = &U::getHash($infoFile);
    my @headers = @{$ihash->{'AffyID'}};
    for (my $i = 0; $i < scalar(@headers); $i++) {
      $ihh->{$headers[$i]} = $i;
    }
  }
  my $vihash;
  my $vihh = {};
  if (-e $vinfoFile) {
    $vihash = &U::getHash($vinfoFile);
    my @headers = @{$vihash->{'ProbeID'}};
    for (my $i = 0; $i < scalar(@headers); $i++) {
      $vihh->{$headers[$i]} = $i;
    }
  }
  
  my @ids = keys(%{$h->{"idHash"}});
  print "AffyID\tName\n";
  for (my $j = 0; $j < scalar(@ids) ; $j++) {
    my $id2 = $ids[$j];
    if (defined $ihash) {
      next if ($ihash->{$id2}->[$ihh->{'hi'}] < $hithr);
      next if ($ihash->{$id2}->[$ihh->{'sd'}] < $sdthr);
      my $dr = $ihash->{$id2}->[$ihh->{'max'}] - $ihash->{$id2}->[$ihh->{'min'}];
      next if ($dr < $drthr);
    }
    if (defined $vihash) {
      my $vdr = $vihash->{$id2}->[$vihh->{'Perc 0.95'}] - $vihash->{$id2}->[$vihh->{'Perc -0.95'}];
      my $sdr1 = $vihash->{$id2}->[$vihh->{'Perc 0.95'}] - $vihash->{$id2}->[$vihh->{'thr'}];
      my $sdr2 = $vihash->{$id2}->[$vihh->{'thr'}] - $vihash->{$id2}->[$vihh->{'Perc -0.95'}];
      my $sdr = abs($sdr1/$sdr2);
      next if (defined $vdrthr && $vdr < $vdrthr);
      next if (defined $sdrthr && $sdrthr > 0 && $sdr < $sdrthr);
      next if (defined $hidrthr && $hidrthr > 0 && $sdr1 < $hidrthr);
      my $sdr = abs($sdr2/$sdr1);
      next if (defined $sdrthr && $sdrthr < 0 && $sdr < -$sdrthr);
      next if (defined $hidrthr && $hidrthr < 0 && $sdr2 < -$hidrthr);
    }
    my $n2 = $h->getName($id2);
    print "$id2\t$n2\n";
  }
}

sub analyzeMerge {
  my (@files) = @_;
  my $keys = [];
  my $hash = {};
  foreach my $f (@files) {
    open(my $fh, "<$f") || die "Can't find $f\n";
    my $num = 0;
    while (<$fh>) {
      s/[\r\n]//g;
      my @list = split("\t");
      if (!defined $hash->{$list[0]}) {
        push @$keys, $list[0];
        push @{$hash->{$list[0]}}, @list;
      }
      else {
        push @{$hash->{$list[0]}}, @list[1 .. $#list];
      }
      if ($num < scalar(@{$hash->{$list[0]}})) {
        $num = scalar(@{$hash->{$list[0]}});
      }
    }
    close($fh);
    foreach my $k (@$keys) {
      if (scalar(@{$hash->{$k}}) < $num) {
        my $n = $num - scalar(@{$hash->{$k}});
        push @{$hash->{$k}}, (map { "" } 1 .. $n);
      }
    }
  }
  foreach my $k (@$keys) {
    print join("\t", @{$hash->{$k}}), "\n";
  }
}

sub analyzeHistogram {
  my ($file, @rest) = @_;
  my $params = {};
  foreach (@rest) {
    my ($k, $v) = split("=");
    $params->{$k} = $v;
  }
  my $res = &U::plotHistogram($file, $params);
  print "min\tmean\tmax\tstddev\tsum\n";
  print join("\t", @{$res->[2]}), "\n";
  print join("\t", "breaks", @{$res->[0]}), "\n";
  print join("\t", "counts", @{$res->[1]}), "\n";
}

sub analyzeHighDR {
  my $pre = shift;
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => $thrf,
      survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $good = [$start .. $end];

  if (1) {
    my @ids = keys(%{$h->{"idHash"}});
    my $infoFile = "$pre\-info.txt";
    my $ihash = &U::getHash($infoFile);
    my @headers = @{$ihash->{'AffyID'}};
    my $ihh = {};
    for (my $i = 0; $i < scalar(@headers); $i++) {
      $ihh->{$headers[$i]} = $i;
    }
    my $count = 0;
    my $data = [];
    foreach my $id (@ids) {
      my $min = $ihash->{$id}->[$ihh->{'min'}];
      my $max = $ihash->{$id}->[$ihh->{'max'}];
      my $dr = $max - $min;
      push @$data, $dr;
    }
    $data = [sort { $a <=> $b } @$data];
    my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
    my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
    print STDERR join("\t", "DR", @$thr), "\n";
    print join("\t", $h->{'headers'}->[0], 
        map { $h->{'headers'}->[$_] } (@$good)), "\n";
    foreach my $id (@ids) {
      my $min = $ihash->{$id}->[$ihh->{'min'}];
      my $max = $ihash->{$id}->[$ihh->{'max'}];
      my $dr = $max - $min;
      next if ($dr < $thr->[0]);
      my $expr = $h->getExprData($id);
      my $e = [ map { if ($_ eq "" || $_ < -2) { -2; } else { $_} }
        map { $expr->[$_] } @$good];
      print join("\t", $id, @$e), "\n";
    }
  }
}

sub analyzeHighDynRange {
  my $pre = shift;
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => "$pre\-survival.txt");
  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $good = [$start .. $end];
  my $f = "$pre\-vinfo.txt";
  open(my $fh, "<$f") || die "Can't open $f\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  my $hhash = {};
  for (my $i = 0; $i < scalar(@headers); $i++) {
    $hhash->{$headers[$i]} = $i;
  }
  my $dhash = {};
  my $data = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $dyn1 = $list[$hhash->{"Perc 0.90"}] - $list[$hhash->{"Perc -0.90"}];
    my $dyn2 = $list[$hhash->{"97.5%"}] - $list[$hhash->{"2.5%"}];
    my $dyn = &U::min([$dyn1, $dyn2]);
    $dhash->{$list[0]} = $dyn;
    push @$data, $dyn;
  }
  close($fh);
  $data = [sort { $a <=> $b } @$data];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  print STDERR join("\t", @$thr), "\n";
  my $ddata = [grep { $_ > $thr->[0] } @$data];
  $ddata = [sort { $a <=> $b } @$ddata];
  my $res = &U::fitStep($ddata, 0, scalar(@$ddata) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  print STDERR join("\t", @$thr), "\n";

  my @ids = keys(%{$h->{"idHash"}});
  print join("\t", $h->{'headers'}->[0], 
      map { $h->{'headers'}->[$_] } (@$good)), "\n";
  foreach my $id (@ids) {
    next if (defined $dhash->{$id} && $dhash->{$id} < $thr->[0]);
    my $expr = $h->getExprData($id);
    my $e = [ map { if ($_ eq "" || $_ < -2) { -2; } else { $_} }
    map { $expr->[$_] } @$good];
    print join("\t", $id, @$e), "\n";
  }
}

sub analyzeHighSDR {
  my $pre = shift;
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => $thrf,
      survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $good = [$start .. $end];

  if (1) {
    my @ids = keys(%{$h->{"idHash"}});
    my $infoFile = "$pre\-info.txt";
    my $ihash = &U::getHash($infoFile);
    my @headers = @{$ihash->{'AffyID'}};
    my $ihh = {};
    for (my $i = 0; $i < scalar(@headers); $i++) {
      $ihh->{$headers[$i]} = $i;
    }
    my $count = 0;
    my $data = [];
    my $sdata = [];
    my $cdata = [];
    foreach my $id (@ids) {
      my $min = $ihash->{$id}->[$ihh->{'min'}];
      my $max = $ihash->{$id}->[$ihh->{'max'}];
      my $dr = $max - $min;
      push @$data, $dr;
      my $sd = $ihash->{$id}->[$ihh->{'sd'}];
      my $hi = $ihash->{$id}->[$ihh->{'hi'}];
      my $lo = $ihash->{$id}->[$ihh->{'lo'}];
      push @$sdata, $sd;
      push @$cdata, &U::min([$hi, $lo]);
    }
    $data = [sort { $a <=> $b } @$data];
    my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
    my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
    print STDERR join("\t", "DR", @$thr), "\n";
    $sdata = [sort { $a <=> $b } @$sdata];
    my $res = &U::fitStep($sdata, 0, scalar(@$sdata) - 1);
    my $sthr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
    print STDERR join("\t", "SD", @$sthr), "\n";
    $cdata = [sort { $a <=> $b } @$cdata];
    my $res = &U::fitStep($cdata, 0, scalar(@$cdata) - 1);
    my $cthr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
    print STDERR join("\t", "C HI LO", @$cthr), "\n";

    print join("\t", $h->{'headers'}->[0], "DR", "SD", "HI", "LO"), "\n";
    foreach my $id (@ids) {
      my $min = $ihash->{$id}->[$ihh->{'min'}];
      my $max = $ihash->{$id}->[$ihh->{'max'}];
      my $sd = $ihash->{$id}->[$ihh->{'sd'}];
      my $hi = $ihash->{$id}->[$ihh->{'hi'}];
      my $lo = $ihash->{$id}->[$ihh->{'lo'}];
      my $dr = $max - $min;
      next if ($dr < $thr->[0]);
      next if ($sd < $sthr->[0]);
      next if ($hi < $cthr->[0]);
      next if ($lo < $cthr->[0]);
      print join("\t", $id, $dr, $sd, $hi, $lo), "\n";
    }
  }
}

sub analyzeCluster {
  my ($file, @rest) = @_;
  my $params = {};
  foreach (@rest) {
    my ($k, $v) = split("=");
    $params->{$k} = $v;
  }
  my $rows = [];
  my $expr = [];
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  my $hhash = {};
  for (my $i = 0; $i < scalar(@headers); $i++) {
    $hhash->{$headers[$i]} = $i;
  }
  my $columns = [ map { $headers[$_] } 1 .. $#headers];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    push @$rows, $list[0];
    push @$expr, [ map { $list[$_] } 1 .. $#headers];
  }
  close($fh);
  my $outfile = "_cluster.png";
  $outfile = $params->{"out"} if (defined $params->{"out"});
  &Hegemon::cluster($outfile, $expr, $rows, $columns, undef, %{$params});
}

sub analyzeBoolean1 {
  my ($file, $num, $st, $p) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    if ($list[$num + 5] > $st && $list[$num + 9] < $p) {
      print join("\t", @list), "\n";
    }
  }
  close($fh);
}

sub convertXLS {
  my $file = shift;
  use Spreadsheet::ParseExcel;
  use Spreadsheet::WriteExcel;
  my $parser = Spreadsheet::ParseExcel->new();
  my $workbook=$parser->parse($file);
  if ( !defined $workbook ) {
    die $parser->error(), ".\n";
  }

  for my $worksheet ( $workbook->worksheets() ) {

    my ( $row_min, $row_max ) = $worksheet->row_range();
    my ( $col_min, $col_max ) = $worksheet->col_range();

    for my $row ( $row_min .. $row_max ) {
      my @l;
      for my $col ( $col_min .. $col_max ) {
        my $cell = $worksheet->get_cell( $row, $col );
        my $v = undef;
        if (defined $cell) {
           $v = $cell->value();
        }
        push @l, $v;
      }
      print join("\t", @l), "\n";
    }
  }
}

sub printProbeInfo {
  my ($file, $gene) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $total = 0;
  my $header = [];
  my $data;
  my $start = 0;
  while (<$fh>) {
    next if (/^#/);
    if (!defined $data && defined $gene && /$gene/i) {
      s/[\r\n]//g;
      my @list = split(/","/);
      $data = [map { s/"//g; $_ } @list];
    }
    if ($start == 0) {
      s/[\r\n]//g;
      my @list = split(/","/);
      $header = [map { s/"//g; $_ } @list];
      $start = 1;
    }
    $total++;
  }
  close($fh);
  print $total, "\n";
  print scalar(@$header), "\n";
  foreach my $i (0 .. $#{$header}) {
    if (defined $data) {
      print "$i: ".$header->[$i].": ".$data->[$i]."\n";
    }
    else {
      print "$i: ".$header->[$i]."\n";
    }
  }
  print "id = 0\n";
  print "name = 0\n";
  print "desc = 0\n";
  print "fields = 0\n";
  print "location = num,chr,start,stop,strand\n";
}

sub analyzeCount {
  my ($type, @dirs) = @_;
  my $expr = [];
  my $arrays = [];
  my $lhash;
  my @queue = (@dirs);
  while (scalar(@queue) > 0) {
    my $dir = shift @queue;
    opendir(my $dirh, $dir) || die "can't opendir $dir: $!";
    while (my $d = readdir($dirh)) {
      next if ($d =~ /^\./);
      if ( -d "$dir/$d" ) {
        push @queue, "$dir/$d";
        next;
      }
      if ($d =~ /^$type$/) {
        my $f = "$dir/$d";
        my @path = split(/\//, $f);
        my $arr1 = $path[$#path - 1];
        if (!defined $lhash) {
          $lhash = {};
          open(my $fh1, "<$f") || die "Can't open $f\n";
          my $h = <$fh1>;
          my $h = <$fh1>; $h =~ s/[\r\n]//g;
          while (<$fh1>) {
            s/[\r\n]//g;
            my @ll = split("\t");
            my $v = $ll[$#ll - 1];
            $lhash->{$ll[0]} = $v;
          }
          close($fh1);
        }
        my $data = {};
        open(my $fh1, "<$f") || die "Can't open $f\n";
        my $h = <$fh1>;
        my $h = <$fh1>; $h =~ s/[\r\n]//g;
        my @hdr = split("\t", $h);
        my $hh1 = {};
        foreach my $i (0 .. $#hdr) { $hh1->{$hdr[$i]} = $i; }
        while (<$fh1>) {
          s/[\r\n]//g;
          my @ll = split("\t");
          my $v = $ll[$#ll];
          $data->{$ll[0]} = $v;
        }
        close($fh1);
        push @$expr, $data;
        push @$arrays, $arr1;
      }
    }
    closedir($dirh);
  }

  my $order = [sort { 
    substr($arrays->[$a], 0, 2) cmp substr($arrays->[$b], 0, 2) ||
      substr($arrays->[$a], 2) <=> substr($arrays->[$b], 2) } 0 .. $#{$arrays}];

  my @ids = keys %{$expr->[0]};
  my $total = [ map { 0 } @$arrays ];
  print join("\t", "ArrayID", "Length", map { $arrays->[$_]} @$order), "\n";
  foreach my $id (@ids) {
    my @e = map { $expr->[$_]->{$id} } @$order;
    $total = [ map { $total->[$_] + $e[$_] } 0 .. $#{$total} ];
    print join("\t", $id, $lhash->{$id}, @e), "\n";
  }
  $total = [ map { sprintf("%.3f", $total->[$_] / 1e6) } 0 .. $#{$total} ];
  print STDERR join("\t", "Total", "-", @$total), "\n";
}

sub analyzeTPM {
  my $f = shift;
  open(my $fh1, "<$f") || die "Can't open $f\n";
  my $h = <$fh1>; $h =~ s/[\r\n]//g;
  my @hdr = split("\t", $h);
  my $hh1 = {};
  my $expr = [];
  my $total = [];
  my $trpk = [];
  foreach my $i (0 .. $#hdr) { $hh1->{$hdr[$i]} = $i; }
  while (<$fh1>) {
    s/[\r\n]//g;
    my @ll = split("\t");
    next if ($ll[1] <= 0);
    my $e = [ map { $ll[$_]/$ll[1] } 2 .. $#ll];
    $trpk = [ map { $e->[$_] + $trpk->[$_] } 0 .. $#{$e}];
    $total = [ map { $ll[$_ + 2] + $total->[$_] } 0 .. $#{$e}];
    push @$expr, [$ll[0], @$e];
  }
  close($fh1);
  print join("\t", $hdr[0], "Name", map { $hdr[$_] } 2 .. $#hdr), "\n";
  foreach my $e (@$expr) {
    print join("\t", $e->[0], $e->[0], map { $e->[$_ + 1]*1e6/$trpk->[$_] } 0 .. $#{$trpk}), "\n";
    #print join("\t", $e->[0], $e->[0], map { $e->[$_ + 1] } 0 .. $#{$trpk}), "\n";
  }
  my $total = [ map { sprintf("%.3f", $total->[$_] / 1e6) } 0 .. $#{$total} ];
  print STDERR join("\t", "Total", "-", @$total), "\n";
}

sub getLog {
  my $v = shift;
  if ($v ne "") {
    if ($v > 1) { $v = log($v)/log(2); }
    else { $v = $v - 1; }
  }
  return $v;
}

sub analyzeTPML {
  my $f = shift;
  open(my $fh1, "<$f") || die "Can't open $f\n";
  my $h = <$fh1>; $h =~ s/[\r\n]//g;
  my @hdr = split("\t", $h);
  my $hh1 = {};
  my $expr = [];
  my $total = [];
  my $trpk = [];
  foreach my $i (0 .. $#hdr) { $hh1->{$hdr[$i]} = $i; }
  while (<$fh1>) {
    s/[\r\n]//g;
    my @ll = split("\t");
    next if ($ll[1] <= 0);
    my $e = [ map { $ll[$_]/$ll[1] } 2 .. $#ll];
    $trpk = [ map { $e->[$_] + $trpk->[$_] } 0 .. $#{$e}];
    $total = [ map { $ll[$_ + 2] + $total->[$_] } 0 .. $#{$e}];
    push @$expr, [$ll[0], @$e];
  }
  close($fh1);
  print join("\t", $hdr[0], "Name", map { $hdr[$_] } 2 .. $#hdr), "\n";
  foreach my $e (@$expr) {
    print join("\t", $e->[0], $e->[0], map { &getLog($_) }  map { $e->[$_ + 1]*1e6/$trpk->[$_] } 0 .. $#{$trpk}), "\n";
    #print join("\t", $e->[0], $e->[0], map { $e->[$_ + 1] } 0 .. $#{$trpk}), "\n";
  }
  my $total = [ map { sprintf("%.3f", $total->[$_] / 1e6) } 0 .. $#{$total} ];
  print STDERR join("\t", "Total", "-", @$total), "\n";
}

sub convertGEO {
  my ($file1, $p, $log, $if, $cf, $ef, $idxf, @rest) = @_;
  my $table1 = &U::getExprSoft($file1, undef, $log);
  my ($platformids, $symbols) = @{$table1->[1]->{$p}};
  my $obj = {};
  $obj->{'hash'} = {};
  $obj->{'headers'} = [];
  $obj->{'cheaders'} = ["Title"];
  $obj->{'platformids'} = $platformids;
  $obj->{'symbols'} = $symbols;
  $obj->{'narray'} = [];
  my $num = scalar(@{$obj->{'cheaders'}});
  $obj->{'carray'} = [0 .. ($num - 1)];

  my $cheader = 0;
  my $chh = {};
  my $hash1 = {};

  foreach my $arr (keys %{$table1->[0]->{$p}}) {
    my $l = $table1->[0]->{$p}->{$arr};
    my $list = [];
    $list->[0] = $l->[0];
    if ($cheader == 0) {
      my @data = map { if (ref($_) ne "HASH") {(split(/: /, $_))[0]} else { ""
} } grep { /: / } @$l;
      my @chdr = map { $data[$_] } 
      grep { length($data[$_]) > 1 && length($data[$_]) <= 80} 
      0 ..  $#data;
      foreach my $ch (@chdr) {
        next if $ch eq "valuefield";
        next if defined $chh->{$ch};
        $chh->{$ch} = scalar(@{$obj->{'cheaders'}});
        push @{$obj->{'cheaders'}}, $ch;
      }
      $num = scalar(@{$obj->{'cheaders'}});
      $obj->{'carray'} = [0 .. ($num - 1)];
      for (my $i = 0; $i < $num; $i++) {
        $chh->{$obj->{'cheaders'}->[$i]} = $i;
      }
    }
  }
  foreach my $arr (keys %{$table1->[0]->{$p}}) {
    my $l = $table1->[0]->{$p}->{$arr};
    my $list = [];
    $list->[0] = $l->[0];
    foreach my $l1 (@$l) {
     foreach my $i (1 .. ($num - 1)) {
       my $str = $obj->{'cheaders'}->[$i];
       $str =~ s/[\(\)\\]/\./g;
       if ($l1 =~ /^$str: (.*)$/) { $list->[$i] = $1; }
     }
    }
    my $title = $l->[0];
    $hash1->{$title} = $arr;
    $list->[$num] = $l->[scalar(@$l) - 1];
    my $time = "";
    my $status = "";
    $obj->{'hash'}->{$arr} = [$time, $status, $list];
    push @{$obj->{'headers'}}, $arr;
  }
  foreach my $i (0 .. ($num - 1)) {
    $obj->{'cheaders'}->[$i] =~ s/\s*\(\s*(.*)\s*\)/-$1/g;
  }

  &U::writeObject($obj, $if, $cf, $ef, $idxf);
}

