#!/usr/bin/perl

if (scalar(@ARGV) < 1) {
    print "Usage: absoluteInfo.pl <cmd> <network.pcl> <network.thr> <network.ph> <phid>\n";
    print "<cmd>: gene pclfile thrfile phfile phid\n";
    print "       list pclfile thrfile listfile phfile phid\n";
    print "       alist <outdir> <arraylist>\n";
    print "       exprPerc pclfile thrfile phfile phid\n";
    print "       dynamicRange pclfile thrfile phfile phid\n";
    print "       dynamicRange-2 pclfile thrfile phfile phid\n";
    print "       dynamicRange-3 pclfile thrfile thr2file phfile phid\n";
    print "       dynamicRange-4 pclfile phfile phid\n";
    print "       dynamicRange-5 pclfile thrfile thr2file phfile phid\n";
    print "       arrayInfo pclfile thrfile arrayid\n";
    print "       listInfo pclfile thrfile listfile\n";
    print "       groupInfo pclfile thrfile listfile\n";
    print "       groupInfo2 pclfile thrfile listfile\n";
    print "       thr pclfile start end gap\n";
    print "       convertpcl exprfile start end\n";
    print "       md5sum listfile igfile\n";
    print "       mergemd5sum listfile1 listfile2 ...\n";
    print "       selectNoDup md5log\n";
    print "       selectPlatform infofile type filelist\n";
    print "       CELinfo listfile\n";
    print "       CELversion listfile\n";
    print "       selectNew listfile1 listfile2 ...\n";
    print "       test\n";
    exit;
}

my $cmd = shift @ARGV;

if ($cmd eq "gene") {
    &printGeneInfo(@ARGV);
}

if ($cmd eq "exprPerc") {
    &printExprPercInfo(@ARGV);
}

if ($cmd eq "dynamicRange") {
    &printDynamicRangeInfo(@ARGV);
}

if ($cmd eq "dynamicRange-2") {
    &printDynamicRangeInfo2(@ARGV);
}

if ($cmd eq "dynamicRange-3") {
    &printDynamicRangeInfo3(@ARGV);
}
if ($cmd eq "dynamicRange-4") {
    &printDynamicRangeInfo4(@ARGV);
}
if ($cmd eq "dynamicRange-5") {
    &printDynamicRangeInfo5(@ARGV);
}
if ($cmd eq "arrayInfo") {
    &printArrayInfo(@ARGV);
}
if ($cmd eq "listInfo") {
    &printListInfo(@ARGV);
}
if ($cmd eq "listInfo1") {
    &printListInfo1(@ARGV);
}
if ($cmd eq "listInfo2") {
    &printListInfo2(@ARGV);
}
if ($cmd eq "groupInfo") {
    &printGroupInfo(@ARGV);
}
if ($cmd eq "groupInfo1") {
    &printGroupInfo1(@ARGV);
}
if ($cmd eq "groupInfo2") {
    &printGroupInfo2(@ARGV);
}
if ($cmd eq "thr") {
    &printThr(@ARGV);
}
if ($cmd eq "convertpcl") {
    &convertpcl(@ARGV);
}
if ($cmd eq "list") {
    &printList(@ARGV);
}
if ($cmd eq "alist") {
    my ($outdir, $listfile, $pclfile, $thrfile, @rest) = @ARGV;
    if (! -e $outdir) {
      mkdir($outdir) || die "Can't create $outdir\n";
    }
    &createGeneInfo($pclfile, $thrfile, $listfile, $outdir);
}

if ($cmd eq "test") {
  my ($min, $max, $mean, $sd, $num, $thrNum) = 
    &getRangeInfo([0, 1, 2, 3, 4, 5.1, 6], 2);
  print "(min, max, mean, sd, num, thrNum)\n";
  print "($min, $max, $mean, $sd, $num, $thrNum)\n";
}
if ($cmd eq "md5sum") {
  &getMd5Sum(@ARGV);
}
if ($cmd eq "mergemd5sum") {
  &mergeMd5Sum(@ARGV);
}
if ($cmd eq "selectNoDup") {
  &selectNoDup(@ARGV);
}
if ($cmd eq "selectPlatform") {
  &selectPlatform(@ARGV);
}
if ($cmd eq "CELinfo") {
  &printCELinfo(@ARGV);
}
if ($cmd eq "CELversion") {
  &printCELversion(@ARGV);
}
if ($cmd eq "selectNew") {
  &selectNew(@ARGV);
}

exit;

sub getPh {
    my ($phfile, $phid) = @_;
    my $res;
    my $fh;
    open($fh, "<$phfile") || die "Can't open $phfile\n";
    while (<$fh>) {
        s/[\r\n]//g;
        my @list = split("\t");
        if ($list[0] eq $phid) {
            $res = [@list];
        }
        if ($phid eq "All") {
            $res = [map {1} @list];
            last;
        }
    }
    close($fh);
    if (!defined $res) {
        print " Undefined Phenotype \n";
        exit;
    }
    return $res;
}

sub printGeneInfo {
    my ($pclfile, $thrfile, $phfile, $phid) = @_;
    my $ph = &getPh($phfile, $phid);
    my $indexHash = {};
    for (my $i = 3; $i < scalar(@$ph); $i++) {
        if ($ph->[$i] == 1) {
            $indexHash->{$i} = 1;
        }
    }
    my $num = scalar(keys(%{$indexHash}));
    if ($num == 0) {
      print STDERR "No matching array for $phid\n"; 
    }
    &printGeneInfoFromIndex(\*STDOUT, $pclfile, $thrfile, $indexHash);
}

sub printGeneInfoFromIndex {
    my ($fh, $pclfile, $thrfile, $indexHash, $idhash) = @_;
    print $fh "AffyID\tname\tthr\tmean\tmean-thr\tperc\tmin\tmax\tsd\n";
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh1>; <$fh1>; <$fh2>;
    while (my $line = <$fh1>) {
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my $id = $list[0];
      $id =~ s/\s//g;
      next if (defined $idhash && !defined $idhash->{$id});
      my $sum = 0;
      my $sum2 = 0;
      my $num = 0;
      my $min = 1000; my $max = 0;
      foreach (keys(%{$indexHash})) {
        $sum += $list[$_];
        $sum2 += $list[$_]*$list[$_];
        if ($min > $list[$_]) {
            $min = $list[$_];
        }
        if ($max < $list[$_]) {
            $max = $list[$_];
        }
        $num ++;
      }
      my $mean = $sum / $num;
      my $mean2 = $sum2 / $num;
      my $sd = sqrt($mean2 - $mean*$mean);
      my $num = 0;
      my $thrNum = 0;
      my $total = scalar(@list);
      for (my $i = 3; $i < scalar(@list); $i++) {
        if ($list[$i] < $mean) {
            $num ++;
        }
        if ($list[$i] < $thr) {
            $thrNum ++;
        }
      }
      if ($num > $thrNum) {
        $total = $total - $thrNum;
      }
      else {
        $total = $thrNum;
      }
      my $perc;
      $perc = ($num - $thrNum) / $total if $total > 0;
      $min =~ s/[\n\r\s]//g;
      $max =~ s/[\n\r\s]//g;
      print $fh $id, "\t", $list[1], "\t", $thr, "\t", $mean, "\t",
            $mean - $thr, "\t", $perc,"\t", $min, "\t", $max, "\t", $sd, "\n";
    }
    close($fh1);
    close($fh2);
}

sub printExprPercInfo {
    my ($pclfile, $thrfile, $phfile, $phid) = @_;
    my $ph = &getPh($phfile, $phid);
    my $indexHash = {};
    for (my $i = 3; $i < scalar(@$ph); $i++) {
        if ($ph->[$i] == 1) {
            $indexHash->{$i} = 1;
        }
    }
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @eweights = split(/\t/, $eweight);
    my @keys = sort keys(%{$indexHash});
    print join("\t", @headers[0..2]), "\t", 
        join("\t", map { $headers[$_] } @keys), "\n";
    print join("\t", @eweights[0..2]), "\t", 
        join("\t", map { $eweights[$_] } @keys), "\n";
    while (my $line = <$fh1>) {
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } @keys;
      my @num = map { 0 } @keys;
      my @thrNum = map { 0 } @keys;
      my @total = map { scalar(@list) } @keys;
      my @perc = map { 0 } @keys;
      for (my $i = 3; $i < scalar(@list); $i++) {
        for (my $j = 0; $j < scalar(@keys); $j++) {
          if ($list[$i] < $expr[$j]) {
            $num[$j] ++;
          }
          if ($list[$i] < $thr) {
            $thrNum[$j] ++;
          }
        }
      }
      for (my $j = 0; $j < scalar(@keys); $j++) {
        if ($num[$j] > $thrNum[$j]) {
          $total[$j] = $total[$j] - $thrNum[$j];
        }
        else {
          $total[$j] = $thrNum[$j];
        }
        $perc[$j] = ($num[$j] - $thrNum[$j]) / $total[$j];
      }
      print join("\t", @list[0..2]), "\t", join("\t", @perc), "\n";
    }
    close($fh1);
    close($fh2);
}

sub getGeneIds {
    my ($listfile) = @_;
    my $res = {};
    my $fh;
    open($fh, "<$listfile") || die "Can't open $listfile\n";
    while (<$fh>) {
        s/[\r\n]//g;
        my @list = split("\t");
        my $id = $list[0];
        $id =~ s/\s//g;
        $res->{$id} = 1;
    }
    close($fh);
    return $res;
}

sub printList {
    my ($pclfile, $thrfile, $listfile, $phfile, $phid) = @_;
    my $idhash = &getGeneIds($listfile);
    my $ph = &getPh($phfile, $phid);
    my $indexHash = {};
    for (my $i = 3; $i < scalar(@$ph); $i++) {
        if ($ph->[$i] == 1) {
            $indexHash->{$i} = 1;
        }
    }
    my $num = scalar(keys(%{$indexHash}));
    if ($num == 0) {
      print STDERR "No matching array for $phid\n"; 
    }
    &printGeneInfoFromIndex(\*STDOUT, $pclfile, $thrfile, $indexHash, $idhash);
}

sub BinarySearch {
   my ($A, $value) = @_;
   my $low = 0;
   my $N = scalar(@$A);
   my $high = $N - 1;
   my $mid = 0;
   while ($low <= $high) {
       $mid = int($low + (($high - $low) / 2));
       if ($A->[$mid] > $value) {
           $high = $mid - 1;
       }
       elsif ($A->[$mid] < $value) {
           $low = $mid + 1;
       }
       else {
           return $mid;
       }
   }
   return $high;
}

sub printDynamicRangeInfo {
    my ($pclfile, $thrfile, $phfile, $phid) = @_;
    my $maxPoints = 100; # -100 to 1, 0, 1, to 100
    my $ph = &getPh($phfile, $phid);
    my $indexHash = {};
    for (my $i = 3; $i < scalar(@$ph); $i++) {
        if ($ph->[$i] == 1) {
            $indexHash->{$i} = 1;
        }
    }
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @eweights = split(/\t/, $eweight);
    my @keys = sort keys(%{$indexHash});
    my @h1 = map { sprintf("Perc -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
    my @h2 = map { sprintf("Perc %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
    my @h3 = map { sprintf("Expr -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
    my @h4 = map { sprintf("Expr %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
    print $headers[0], "\t", join("\t", @h1), "\t", join("\t", @h2);
    print              "\t", join("\t", @h3), "\t", join("\t", @h4), "\n";
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } @keys;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      # Debug prints
      #print scalar(@sortedExprBelowThr), "\n";
      #print scalar(@sortedExprAboveThr), "\n";
      #print $thr, "\n";
      #print &BinarySearch([ 1, 2, 3, 4, 4, 5, 5, 6, 8, 9, 10], 5), "\n";
      print $list[0];
      for (my $i = 0; $i < $maxPoints; $i++) {
        if (scalar(@sortedExprBelowThr) != 0) {
          my $index = int($i*1.0*$#sortedExprBelowThr/$maxPoints);
          my $expr = $sortedExprBelowThr[$index];
          print "\t", sprintf("%.2f", $expr);
        }
        else {
          print "\t", 0;
        }
      }
      for (my $i = 0; $i <= $maxPoints; $i++) {
        if (scalar(@sortedExprAboveThr) != 0) {
          my $index = int($i*1.0*$#sortedExprAboveThr/$maxPoints);
          my $expr = $sortedExprAboveThr[$index];
          if ($i == 0) {
            $expr = $thr;
          }
          print "\t", sprintf("%.2f", $expr);
        }
        else {
          print "\t", 0;
        }
      }
      my $min_below = $sortedExprBelowThr[0];
      my $max_below = $sortedExprBelowThr[$#sortedExprBelowThr];
      for (my $i = 0; $i < $maxPoints; $i++) {
        if (scalar(@sortedExprBelowThr) != 0) {
          my $expr = $min_below + $i*($max_below - $min_below)/$maxPoints;
          my $index = &BinarySearch(\@sortedExprBelowThr, $expr);
          my $num = scalar(@sortedExprBelowThr);
          my $perc = ($num - $index)*1.0/$num;
          print "\t", sprintf("-%.2f", $perc);
        }
        else {
          print "\t", 1;
        }
      }
      my $min_above = $sortedExprAboveThr[0];
      my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
      for (my $i = 0; $i <= $maxPoints; $i++) {
        if (scalar(@sortedExprAboveThr) != 0) {
          my $expr = $min_above + $i*($max_above - $min_above)/$maxPoints;
          my $index = &BinarySearch(\@sortedExprAboveThr, $expr);
          my $num = scalar(@sortedExprAboveThr);
          my $perc = $index*1.0/$num;
          print "\t", sprintf("%.2f", $perc);
        }
        else {
          print "\t", 1;
        }
      }
      print "\n";
    }
    close($fh1);
    close($fh2);
}

sub printDynamicRangeInfo2 {
    my ($pclfile, $thrfile, $phfile, $phid) = @_;
    my $maxPoints = 100; # -100 to 1, 0, 1, to 100
    my $ph = &getPh($phfile, $phid);
    my $indexHash = {};
    for (my $i = 3; $i < scalar(@$ph); $i++) {
        if ($ph->[$i] == 1) {
            $indexHash->{$i} = 1;
        }
    }
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @eweights = split(/\t/, $eweight);
    my @keys = sort keys(%{$indexHash});
    my @h1 = map { sprintf("Perc -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
    my @h2 = map { sprintf("Perc %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
    my @h3 = map { sprintf("Expr -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
    my @h4 = map { sprintf("Expr %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
    print $headers[0], "\t", join("\t", @h1), "\t", join("\t", @h2);
    print              "\t", join("\t", @h3), "\t", join("\t", @h4), "\n";
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } @keys;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      # Debug prints
      #print scalar(@sortedExprBelowThr), "\n";
      #print scalar(@sortedExprAboveThr), "\n";
      #print $thr, "\n";
      #print &BinarySearch([ 1, 2, 3, 4, 4, 5, 5, 6, 8, 9, 10], 5), "\n";
      print $list[0];
      for (my $i = 0; $i < $maxPoints; $i++) {
        if (scalar(@sortedExprBelowThr) != 0) {
          my $index = int($i*1.0*$#sortedExprBelowThr/$maxPoints);
          my $expr = $sortedExprBelowThr[$index];
          print "\t", sprintf("%.2f", $expr);
        }
        else {
          print "\t", 0;
        }
      }
      for (my $i = 0; $i <= $maxPoints; $i++) {
        if (scalar(@sortedExprAboveThr) != 0) {
          my $index = int($i*1.0*$#sortedExprAboveThr/$maxPoints);
          my $expr = $sortedExprAboveThr[$index];
          if ($i == 0) {
            $expr = $thr;
          }
          print "\t", sprintf("%.2f", $expr);
        }
        else {
          print "\t", 0;
        }
      }
      my $min_below = $sortedExprBelowThr[0];
      my $max_below = $sortedExprBelowThr[$#sortedExprBelowThr];
      for (my $i = 0; $i < $maxPoints; $i++) {
        if (scalar(@sortedExprBelowThr) != 0) {
          my $expr = $min_below + $i*($max_below - $min_below)/$maxPoints;
          my $index = &BinarySearch(\@sortedExprBelowThr, $expr);
          my $num = scalar(@sortedExprBelowThr);
          print "\t", sprintf("%d", $index);
        }
        else {
          print "\t", 1;
        }
      }
      my $min_above = $sortedExprAboveThr[0];
      my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
      for (my $i = 0; $i <= $maxPoints; $i++) {
        if (scalar(@sortedExprAboveThr) != 0) {
          my $expr = $min_above + $i*($max_above - $min_above)/$maxPoints;
          my $index = &BinarySearch(\@sortedExprAboveThr, $expr);
          my $num = scalar(@sortedExprAboveThr);
          my $perc = $index + scalar(@sortedExprBelowThr);
          print "\t", sprintf("%d", $perc);
        }
        else {
          print "\t", 1;
        }
      }
      print "\n";
    }
    close($fh1);
    close($fh2);
}

sub getHeaders {
    my $pclfile = shift;
    my $fh1;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    my $header=<$fh1>; 
    $header =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    close($fh1);
    return [@headers];
}

sub createGeneInfo {
    my ($pclfile, $thrfile, $listfile, $outdir) = @_;
    my $headers = &getHeaders($pclfile);
    my $hHash = {};
    for (my $i = 0; $i < scalar(@$headers); $i++) {
      $hHash->{$headers->[$i]} = $i;
    }
    my $fh;
    open($fh, "<$listfile") || die "Can't open $listfile";
    while (<$fh>) {
      s/[\r\n]//g;
      my @list = split("\t");
      my $indexHash = {};
      foreach (@list) {
        if (defined $hHash->{$_}) {
          $indexHash->{$hHash->{$_}} = 1;
        }
      }
      my $name = lc($list[0]);
      $name =~ s/[-\s\+\(\)]/_/g;
      my $num = scalar(keys(%{$indexHash}));
      print STDERR "$num matching array for $name\n"; 
      if ($num == 0) {
        next;
      }
      my $ofh;
      open($ofh, ">$outdir/$name\_geneinfo.txt") || die "Can't open $outdir/$name\_geneinfo.txt\n";
      &printGeneInfoFromIndex($ofh, $pclfile, $thrfile, $indexHash);
      close($ofh);
    }
    close($fh);
}

sub getRangeInfo {
  my ($expr, $thr) = @_;
  my $sum = 0;
  my $sum2 = 0;
  my $num = 0;
  my $thrNum = 0;
  my $min = $expr->[0]; my $max = $expr->[0];
  foreach (@$expr) {
    next if (/^\s*$/);
    $sum += $_;
    $sum2 += $_*$_;
    if ($min > $_) {
      $min = $_;
    }
    if ($max < $_) {
      $max = $_;
    }
    if ($_ < $thr) {
      $thrNum ++;
    }
    $num ++;
  }
  my $mean = $sum / $num;
  my $mean2 = $sum2 / $num;
  my $sd = sqrt($mean2 - $mean*$mean);
  return ($min, $max, $mean, $sd, $num, $thrNum);
}

sub printDynamicRangeInfo3 {
    my ($pclfile, $thrfile, $thr2file, $phfile, $phid) = @_;
    my $maxPoints = 20; # -100 to 1, 0, 1, to 100
    my $gmin = 1.0;
    my $gmax = 16.0;
    my $gstep = 0.2;
    my $gnum = int(($gmax - $gmin)/$gstep);
    my $ph = &getPh($phfile, $phid);
    my $indexHash = {};
    for (my $i = 3; $i < scalar(@$ph); $i++) {
        if ($ph->[$i] == 1) {
            $indexHash->{$i} = 1;
        }
    }
    my $fh1;
    my $fh2;
    my $fh3;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    open($fh3, "<$thr2file") || die "Can't open $thr2file";
    <$fh2>; <$fh3>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @eweights = split(/\t/, $eweight);
    my @keys = sort keys(%{$indexHash});
    my @h1 = map { sprintf("Perc -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
    my @h2 = map { sprintf("Perc %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
    my @h3 = map { sprintf("Expr -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
    my @h4 = map { sprintf("Expr %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
    my @h5 = map { sprintf("Abs %.2f", $_*$gstep + $gmin) } 0 .. $gnum;
    print join("\t", $headers[0], "min", "max", "mean", "sd", "thr", "thrNum", "pval");
    print "\t", join("\t", @h1), "\t", join("\t", @h2);
    print "\t", join("\t", @h3), "\t", join("\t", @h4);
    print "\t", join("\t", @h5), "\n";
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my $thr2Line = <$fh3>; $thr2Line =~ s/[\r\n]//g;
      my ($n, $pval, $rest) = split(/\t/, $thr2Line);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } @keys;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my ($min, $max, $mean, $sd, $num, $thrNum) = &getRangeInfo(\@expr, $thr);
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      # Debug prints
      #print scalar(@sortedExprBelowThr), "\n";
      #print scalar(@sortedExprAboveThr), "\n";
      #print $thr, "\n";
      #print &BinarySearch([ 1, 2, 3, 4, 4, 5, 5, 6, 8, 9, 10], 5), "\n";
      print join("\t", $list[0], $min, $max, $mean, $sd, $thr, $num - $thrNum, $pval);
      for (my $i = 0; $i < $maxPoints; $i++) {
        if (scalar(@sortedExprBelowThr) != 0) {
          my $index = int($i*1.0*$#sortedExprBelowThr/$maxPoints);
          my $expr = $sortedExprBelowThr[$index];
          print "\t", sprintf("%.2f", $expr);
        }
        else {
          print "\t", 0;
        }
      }
      for (my $i = 0; $i <= $maxPoints; $i++) {
        if (scalar(@sortedExprAboveThr) != 0) {
          my $index = int($i*1.0*$#sortedExprAboveThr/$maxPoints);
          my $expr = $sortedExprAboveThr[$index];
          if ($i == 0) {
            $expr = $thr;
          }
          print "\t", sprintf("%.2f", $expr);
        }
        else {
          print "\t", 0;
        }
      }
      my $min_below = $sortedExprBelowThr[0];
      my $max_below = $sortedExprBelowThr[$#sortedExprBelowThr];
      my $last_index = 0;
      for (my $i = 0; $i < $maxPoints; $i++) {
        if (scalar(@sortedExprBelowThr) != 0) {
          my $expr = $min_below + $i*($max_below - $min_below)/$maxPoints;
          my $index = &BinarySearch(\@sortedExprBelowThr, $expr);
          print "\t", sprintf("%d", $index - $last_index);
          $last_index = $index;
        }
        else {
          print "\t", 1;
        }
      }
      my $min_above = $sortedExprAboveThr[0];
      my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
      for (my $i = 0; $i <= $maxPoints; $i++) {
        if (scalar(@sortedExprAboveThr) != 0) {
          my $expr = $min_above + $i*($max_above - $min_above)/$maxPoints;
          my $index = &BinarySearch(\@sortedExprAboveThr, $expr);
          my $new_index = $index + scalar(@sortedExprBelowThr);
          print "\t", sprintf("%d", $new_index - $last_index);
          $last_index = $new_index;
        }
        else {
          print "\t", 1;
        }
      }
      my $last_index = 0;
      for (my $i = 0; $i <= $gnum; $i++) {
        if (scalar(@sortedExpr) != 0) {
          my $expr = $gmin + $i*$gstep;
          my $index = &BinarySearch(\@sortedExpr, $expr);
          if ($index < 0) {
            $index = 0;
          }
          print "\t", sprintf("%d", $index - $last_index);
          $last_index = $index;
        }
        else {
          print "\t", 0;
        }
      }
      print "\n";
    }
    close($fh1);
    close($fh2);
    close($fh3);
}

sub printArrayInfo {
    my ($pclfile, $thrfile, $aid) = @_;
    my $maxPoints = 20; # -100 to 1, 0, 1, to 100
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my $hhash = {};
    for (my $i = 3; $i < scalar(@headers); $i++) {
      $hhash->{$headers[$i]} = $i;
    }
    if (!defined $hhash->{$aid}) {
      print STDERR "$aid not found\n";
      return;
    }
    my $lineno = 0;
    my $pos = $hhash->{$aid};
    print "AffyID\tExpr\tPerc\tRatio\n";
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } 3 .. $#headers;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      my $exp = $list[$pos];
      my ($perc, $ratio);
      if ($exp < $thr) {
        my $index = &BinarySearch(\@sortedExprBelowThr, $exp);
        my $min_below = $sortedExprBelowThr[0];
        $perc = ($index - $#sortedExprBelowThr) * 1.0 /$#sortedExprBelowThr;
        $ratio = ($exp - $thr) / ($thr - $min_below);
      }
      else {
        my $index = &BinarySearch(\@sortedExprAboveThr, $exp);
        my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
        $perc = $index * 1.0 /$#sortedExprAboveThr;
        $ratio = ($exp - $thr) / ($max_above - $thr);
      }
      print join("\t", $list[0], $exp, $perc, $ratio), "\n";
    }
    close($fh1);
    close($fh2);
}

sub printListInfo {
    my ($pclfile, $thrfile, $listfile) = @_;
    my $fh;
    open($fh, "<$listfile") || die "Can't open $listfile";
    my $idHash = {};
    while (<$fh>) {
      s/[\r\n]//g;
      my @list = split(/\//);
      my $name = $list[$#list];
      $name =~ s/.CEL$//ig;
      $idHash->{$name} = 1;
    }
    close($fh);

    my $maxPoints = 20; # -100 to 1, 0, 1, to 100
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my $found = 0;
    for (my $i = 3; $i < scalar(@headers); $i++) {
      if (defined $idHash->{$headers[$i]}) {
        $idHash->{$headers[$i]} = $i;
        $found ++;
      }
    }
    if ($found <= 0) {
      print STDERR "arrays not found\n";
      return;
    }
    my @fhs;
    my @names = keys %{$idHash};
    my $i = 0;
    foreach my $f (@names) {
      open($fhs[$i], ">info/$f\-info.txt") || die "Can't open $f\-info.txt\n";
      my $fh = $fhs[$i];
      print $fh "AffyID\tExpr\tPerc\tRatio\n";
      $i++;
    }
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } 3 .. $#headers;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      my $ii = 0;
      foreach $f (@names) {
        my $pos = $idHash->{$f};
        my $exp = $list[$pos];
        my ($perc, $ratio);
        if ($exp < $thr) {
          my $index = &BinarySearch(\@sortedExprBelowThr, $exp);
          my $min_below = $sortedExprBelowThr[0];
          $perc = ($index - $#sortedExprBelowThr) * 1.0 /$#sortedExprBelowThr;
          $ratio = ($exp - $thr) / ($thr - $min_below);
        }
        else {
          my $index = &BinarySearch(\@sortedExprAboveThr, $exp);
          my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
          $perc = $index * 1.0 /$#sortedExprAboveThr;
          $ratio = ($exp - $thr) / ($max_above - $thr);
        }
        my $fh = $fhs[$ii++];
        print $fh join("\t", $list[0], $exp, $perc, $ratio), "\n";
      }
    }
    close($fh1);
    close($fh2);
    my $i = 0;
    foreach my $f (@names) {
      close($fhs[$i++]);
    }
}

sub printListInfo1 {
    my ($pclfile, $thrfile, $listfile) = @_;
    my $fh;
    open($fh, "<$listfile") || die "Can't open $listfile";
    my $groupHash = {};
    my $idHash = {};
    while (<$fh>) {
      $f = $_;
      $f =~ s/[\r\n]//g;
      my @list = split(/\//, $f);
      my $name = $list[$#list];
      $name =~ s/.pcl$//ig;
      open(my $fp, "<$f") || die "Can't open $f\n";
      my $head = <$fp>;
      my $ew = <$fp>;
      $head =~ s/[\r\n]//g;
      my @headers = split("\t", $head);
      $groupHash->{$name} = [$fp, [@headers]];
      for (my $i = 3; $i < scalar(@headers); $i++) {
        $idHash->{$headers[$i]} = [$name, $fp, $i];
      }
    }
    close($fh);
    if (scalar(keys %{$idHash}) <= 0) {
      print STDERR "arrays not found\n";
      return;
    }

    my $maxPoints = 20; # -100 to 1, 0, 1, to 100
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @fhs;
    my @names = keys %{$idHash};
    my $i = 0;
    foreach my $f (@names) {
      open($fhs[$i], ">info/$f\-info.txt") || die "Can't open $f\-info.txt\n";
      my $fh = $fhs[$i];
      print $fh "AffyID\tExpr\tPerc\tRatio\n";
      $i++;
    }
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } 3 .. $#headers;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      my $gData = {};
      foreach my $g (keys %{$groupHash}) {
        my $fp = $groupHash->{$g}->[0];
        my $line = <$fp>;
        $line =~ s/[\r\n]//g;
        my @l = split("\t", $line);
        $gData->{$g} = [@l];
      }
      my $ii = 0;
      foreach $f (@names) {
        my $pos = $idHash->{$f}->[2];
        my $exp = $gData->{$idHash->{$f}->[0]}->[$pos];
        my ($perc, $ratio);
        if ($exp < $thr) {
          my $index = &BinarySearch(\@sortedExprBelowThr, $exp);
          my $min_below = $sortedExprBelowThr[0];
          $perc = ($index - $#sortedExprBelowThr) * 1.0 /$#sortedExprBelowThr;
          $ratio = ($exp - $thr) / ($thr - $min_below);
        }
        else {
          my $index = &BinarySearch(\@sortedExprAboveThr, $exp);
          my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
          $perc = $index * 1.0 /$#sortedExprAboveThr;
          $ratio = ($exp - $thr) / ($max_above - $thr);
        }
        my $fh = $fhs[$ii++];
        print $fh join("\t", $list[0], $exp, $perc, $ratio), "\n";
      }
    }
    close($fh1);
    close($fh2);
    my $i = 0;
    foreach my $f (@names) {
      close($fhs[$i++]);
    }
    foreach my $g (keys %{$groupHash}) {
      my $fp = $groupHash->{$g}->[0];
      close($fp);
    }
}

sub printListInfo2 {
    my ($pclfile, $thrfile, $listfile, @files) = @_;
    my $fhHash = {};
    my @fileps;
    for (my $i = 0; $i < scalar(@files); $i++) {
      my $f = $files[$i];
      open(my $fp, "<$f") || die "Can't open $f\n";
      $fileps[$i] = $fp;
      my $head = <$fp>;
      my $ew = <$fp>;
      $head =~ s/[\r\n]//g;
      my @headers = split("\t", $head);
      for (my $j = 1; $j < scalar(@headers); $j++) {
        $fhHash->{$headers[$j]} = [$fp, \@headers, $i, $j];
      }
    }
    my $fh;
    open($fh, "<$listfile") || die "Can't open $listfile";
    my $groupHash = {};
    my $idHash = {};
    while (<$fh>) {
      s/[\r\n]//g;
      my ($gname, $h) = split("\t", $_);
      if (!defined $fhHash->{$h}) {
        print STDERR "$h not found\n";
        exit(1);
      }
      my @list = split("_", $gname);
      pop @list;
      my $name = join("_", @list);
      if (!defined $groupHash->{$name}) {
        $groupHash->{$name} = [$fhHash->{$h}];
      }
      else {
        push @{$groupHash->{$name}}, $fhHash->{$h};
      }
      if (defined $idHash->{$gname}) {
        print STDERR "$gname found in duplicate\n";
        exit(1);
      }
      $idHash->{$gname} = [$name, @{$fhHash->{$h}}];
    }
    close($fh);
    if (scalar(keys %{$idHash}) <= 0) {
      print STDERR "arrays not found\n";
      return;
    }

    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @fhs;
    my @names = keys %{$idHash};
    my $i = 0;
    foreach my $f (@names) {
      open($fhs[$i], ">info/$f\-info.txt") || die "Can't open $f\-info.txt\n";
      my $fh = $fhs[$i];
      print $fh "AffyID\tExpr\tPerc\tRatio\n";
      $i++;
    }
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } 3 .. $#headers;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      my $gData = {};
      for (my $i = 0; $i < scalar(@fileps); $i++) {
        my $fp = $fileps[$i];
        my $line = <$fp>;
        $line =~ s/[\r\n]//g;
        my @l = split("\t", $line);
        $gData->{$i} = [@l];
      }
      my $ii = 0;
      foreach $f (@names) {
        my $pos = $idHash->{$f}->[4];
        my $exp = $gData->{$idHash->{$f}->[3]}->[$pos];
        my ($perc, $ratio);
        if ($exp < $thr) {
          my $index = &BinarySearch(\@sortedExprBelowThr, $exp);
          my $min_below = $sortedExprBelowThr[0];
          $perc = ($index - $#sortedExprBelowThr) * 1.0 /$#sortedExprBelowThr;
          $ratio = ($exp - $thr) / ($thr - $min_below);
        }
        else {
          my $index = &BinarySearch(\@sortedExprAboveThr, $exp);
          my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
          $perc = $index * 1.0 /$#sortedExprAboveThr;
          $ratio = ($exp - $thr) / ($max_above - $thr);
        }
        my $fh = $fhs[$ii++];
        print $fh join("\t", $list[0], $exp, $perc, $ratio), "\n";
      }
    }
    close($fh1);
    close($fh2);
    my $i = 0;
    foreach my $f (@names) {
      close($fhs[$i++]);
    }
    for (my $i = 0; $i < scalar(@fileps); $i++) {
      my $fp = $fileps[$i];
      close($fp);
    }
}

sub printGroupInfo1 {
    my ($pclfile, $thrfile, $listfile) = @_;
    my $fh;
    open($fh, "<$listfile") || die "Can't open $listfile";
    my $groupHash = {};
    my $idHash = {};
    while (<$fh>) {
      $f = $_;
      $f =~ s/[\r\n]//g;
      my @list = split(/\//, $f);
      my $name = $list[$#list];
      $name =~ s/.pcl$//ig;
      open(my $fp, "<$f") || die "Can't open $f\n";
      my $head = <$fp>;
      my $ew = <$fp>;
      $head =~ s/[\r\n]//g;
      my @headers = split("\t", $head);
      $groupHash->{$name} = [$fp, [@headers]];
      for (my $i = 3; $i < scalar(@headers); $i++) {
        $idHash->{$headers[$i]} = [$name, $fp, $i];
      }
    }
    close($fh);
    if (scalar(keys %{$idHash}) <= 0) {
      print STDERR "arrays not found\n";
      return;
    }

    my $maxPoints = 20; # -100 to 1, 0, 1, to 100
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @fhs;
    my @names = keys %{$groupHash};
    my $i = 0;
    foreach my $f (@names) {
      open($fhs[$i], ">info/$f\-info.txt") || die "Can't open $f\-info.txt\n";
      my $fh = $fhs[$i];
      print $fh "AffyID\tExpr\tPerc\tRatio\n";
      $i++;
    }
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } 3 .. $#headers;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      my $ii = 0;
      foreach $g (@names) {
        my $fp = $groupHash->{$g}->[0];
        my $line = <$fp>;
        $line =~ s/[\r\n]//g;
        my @l = split("\t", $line);
        my @exprs = @l[3 .. $#l];
        my $exp = &mean(\@exprs, 0, $#exprs);
        my ($perc, $ratio);
        if ($exp < $thr) {
          my $index = &BinarySearch(\@sortedExprBelowThr, $exp);
          my $min_below = $sortedExprBelowThr[0];
          $perc = ($index - $#sortedExprBelowThr) * 1.0 /$#sortedExprBelowThr;
          $ratio = ($exp - $thr) / ($thr - $min_below);
        }
        else {
          my $index = &BinarySearch(\@sortedExprAboveThr, $exp);
          my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
          $perc = $index * 1.0 /$#sortedExprAboveThr;
          $ratio = ($exp - $thr) / ($max_above - $thr);
        }
        my $fh = $fhs[$ii++];
        print $fh join("\t", $list[0], $exp, $perc, $ratio), "\n";
      }
    }
    close($fh1);
    close($fh2);
    my $i = 0;
    foreach my $f (@names) {
      close($fhs[$i++]);
    }
    foreach my $g (keys %{$groupHash}) {
      my $fp = $groupHash->{$g}->[0];
      close($fp);
    }
}

sub printGroupInfo {
    my ($pclfile, $thrfile, $listfile) = @_;
    my $fh;
    open($fh, "<$listfile") || die "Can't open $listfile";
    my $idHash = {};
    while (<$fh>) {
      s/[\r\n]//g;
      my @list = split(/\//);
      my $name = $list[$#list];
      $name =~ s/.CEL$//ig;
      $idHash->{$name} = 1;
    }
    close($fh);

    my $maxPoints = 20; # -100 to 1, 0, 1, to 100
    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my $found = 0;
    for (my $i = 3; $i < scalar(@headers); $i++) {
      if (defined $idHash->{$headers[$i]}) {
        $idHash->{$headers[$i]} = $i;
        $found ++;
      }
    }
    if ($found <= 0) {
      print STDERR "arrays not found\n";
      return;
    }
    my @fhs;
    my @names = keys %{$idHash};
    my $groupHash = {};
    foreach my $f (@names) {
       my @list = split("_", $f); 
       my $end = $#list - 1;
       my $g = join("_", @list[0 .. $end]);
       push @{$groupHash->{$g}}, $f;
    }
    my @groups = keys %{$groupHash};

    my $i = 0;
    foreach my $f (@groups) {
      open($fhs[$i], ">info/$f\-info.txt") || die "Can't open $f\-info.txt\n";
      my $fh = $fhs[$i];
      print $fh "AffyID\tExpr\tPerc\tRatio\n";
      $i++;
    }
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } 3 .. $#headers;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      my $ii = 0;
      foreach $g (@groups) {
        my @exprs = map { $list[$idHash->{$_}] } @{$groupHash->{$g}};
        my $exp = &mean(\@exprs, 0, $#exprs);
        my ($perc, $ratio);
        if ($exp < $thr) {
          my $index = &BinarySearch(\@sortedExprBelowThr, $exp);
          my $min_below = $sortedExprBelowThr[0];
          $perc = ($index - $#sortedExprBelowThr) * 1.0 /$#sortedExprBelowThr;
          $ratio = ($exp - $thr) / ($thr - $min_below);
        }
        else {
          my $index = &BinarySearch(\@sortedExprAboveThr, $exp);
          my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
          $perc = $index * 1.0 /$#sortedExprAboveThr;
          $ratio = ($exp - $thr) / ($max_above - $thr);
        }
        my $fh = $fhs[$ii++];
        print $fh join("\t", $list[0], $exp, $perc, $ratio), "\n";
      }
    }
    close($fh1);
    close($fh2);
    my $i = 0;
    foreach my $f (@groups) {
      close($fhs[$i++]);
    }
}

sub printGroupInfo2 {
    my ($pclfile, $thrfile, $listfile, @files) = @_;
    my $fhHash = {};
    my @fileps;
    for (my $i = 0; $i < scalar(@files); $i++) {
      my $f = $files[$i];
      open(my $fp, "<$f") || die "Can't open $f\n";
      $fileps[$i] = $fp;
      my $head = <$fp>;
      my $ew = <$fp>;
      $head =~ s/[\r\n]//g;
      my @headers = split("\t", $head);
      for (my $j = 1; $j < scalar(@headers); $j++) {
        $fhHash->{$headers[$j]} = [$fp, \@headers, $i, $j];
      }
    }
    my $fh;
    open($fh, "<$listfile") || die "Can't open $listfile";
    my $groupHash = {};
    my $idHash = {};
    while (<$fh>) {
      s/[\r\n]//g;
      my ($gname, $h) = split("\t", $_);
      if (!defined $fhHash->{$h}) {
        print STDERR "$h not found\n";
        exit(1);
      }
      my @list = split("_", $gname);
      pop @list;
      my $name = join("_", @list);
      if (!defined $groupHash->{$name}) {
        $groupHash->{$name} = [$fhHash->{$h}];
      }
      else {
        push @{$groupHash->{$name}}, $fhHash->{$h};
      }
      if (defined $idHash->{$gname}) {
        print STDERR "$gname found in duplicate\n";
        exit(1);
      }
      $idHash->{$gname} = [$name, @{$fhHash->{$h}}];
    }
    close($fh);
    if (scalar(keys %{$idHash}) <= 0) {
      print STDERR "arrays not found\n";
      return;
    }

    my $fh1;
    my $fh2;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    <$fh2>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @fhs;
    my @names = keys %{$groupHash};
    my $i = 0;
    foreach my $f (@names) {
      open($fhs[$i], ">info/$f\-info.txt") || die "Can't open $f\-info.txt\n";
      my $fh = $fhs[$i];
      print $fh "AffyID\tExpr\tPerc\tRatio\n";
      $i++;
    }
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } 3 .. $#headers;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my $min = $sortedExpr[0];
      my $max = $sortedExpr[$#sortedExpr];
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      my $gData = {};
      for (my $i = 0; $i < scalar(@fileps); $i++) {
        my $fp = $fileps[$i];
        my $line = <$fp>;
        $line =~ s/[\r\n]//g;
        my @l = split("\t", $line);
        $gData->{$i} = [@l];
      }
      my $ii = 0;
      foreach $g (@names) {
        my @exprs = map { $gData->{$_->[2]}->[$_->[3]] } @{$groupHash->{$g}};
        my $exp = &mean(\@exprs, 0, $#exprs);
        my ($perc, $ratio);
        if ($exp < $thr) {
          my $index = &BinarySearch(\@sortedExprBelowThr, $exp);
          my $min_below = $sortedExprBelowThr[0];
          $perc = ($index - $#sortedExprBelowThr) * 1.0 /$#sortedExprBelowThr;
          $ratio = ($exp - $thr) / ($thr - $min_below);
        }
        else {
          my $index = &BinarySearch(\@sortedExprAboveThr, $exp);
          my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
          $perc = $index * 1.0 /$#sortedExprAboveThr;
          $ratio = ($exp - $thr) / ($max_above - $thr);
        }
        my $fh = $fhs[$ii++];
        print $fh join("\t", $list[0], $exp, $perc, $ratio), "\n";
      }
    }
    close($fh1);
    close($fh2);
    my $i = 0;
    foreach my $f (@names) {
      close($fhs[$i++]);
    }
    for (my $i = 0; $i < scalar(@fileps); $i++) {
      my $fp = $fileps[$i];
      close($fp);
    }
}

sub mse {
  my ($arrayref, $start, $end) = @_;
  my $m = &mean($arrayref, $start, $end);
  my $result = 0.0;
  for (my $i = $start; $i <= $end; $i++) {
    $result += ($arrayref->[$i] - $m) * ($arrayref->[$i] - $m);
  }
  return $result;
}

sub mean {
  my ($arrayref, $start, $end) = @_;
  my $result = 0.0;
  if ($start > $end) {
    return $result;
  }
  for (my $i = $start; $i <= $end; $i++) {
    $result += $arrayref->[$i];
  }
  return $result / ($end - $start + 1);
}

sub sum {
  my ($arrayref, $start, $end) = @_;
  my $result = 0.0;
  for (my $i = $start; $i <= $end; $i++) {
    $result += $arrayref->[$i];
  }
  return $result;
}

sub fitStep {
  my ($data, $start, $end) = @_;
  my $count = $end - $start + 1;
  if ($count <= 0) {
    return [0, 0, 0, 0, 0, 0, 0, 0];
  }
  my @sseArray = map { 0.0 } 0 .. ($count - 1);
  my $sum = &sum($data, $start, $end);
  my $mean = &mean($data, $start, $end);
  my $sstot = &mse($data, $start, $end);
  my $sum1 = 0.0;
  my $count1 = 0;
  my $m1 = 0.0;
  my $sum2 = $sum;
  my $count2 = $count;
  my $m2 = ($sum/$count);
  my $sum1sq = 0.0;
  my $sum2sq = $sstot;
  my $sse = $sum1sq + $sum2sq;
  for (my $i = 0; $i < $count; $i++) {
    my $entry = $data->[$i + $start];
    if (!defined $entry || $entry eq "") {
      $sseArray[$i] = $sse;
      next;
    }
    $count1 ++;
    $count2 --;
    if ($count2 == 0) {
      $sseArray[$i] = $sstot;
      next;
    }
    my $tmp = ($mean - ($entry + $sum1)/$count1);
    $sum1sq = $sum1sq + ($entry - $mean) * ($entry - $mean) - 
      $tmp * $tmp * $count1 + ($count1 - 1) * ($mean - $m1) * ($mean - $m1);
    $tmp = ($mean - ($sum2 - $entry)/$count2);
    $sum2sq = $sum2sq - ($entry - $mean) * ($entry - $mean) - 
      $tmp * $tmp * $count2 + ($count2 + 1) * ($mean - $m2) * ($mean - $m2);
    $sum1 += $entry;
    $sum2 -= $entry;
    $m1 = $sum1/$count1;
    $m2 = $sum2/$count2;
    $sse = $sum1sq + $sum2sq;
    $sseArray[$i] = $sse;
  }
  #print join(",", map { sprintf("%.2f", $_) } @sseArray), "\n-----------\n";
  my $bestSse;
  my $bestIndex = 0;
  for (my $i = 0; $i < $count ; $i++) {
    my $index = $i + $start;
    if (!defined $bestSse) {
      $bestSse = $sseArray[$i];
      $bestIndex = $index;
    }
    if ($sseArray[$i] < $bestSse) {
      $bestSse = $sseArray[$i];
      $bestIndex = $index;
    }
    #{# Debug code start
    #  $entry = $data->[$index];
    #  $sum1sq = &mse($data, $start, $index);
    #  $sum2sq = &mse($data, $index + 1, $end);
    #  $sse = $sum1sq + $sum2sq;
    #  printf "%.2f\t%.2f\t%.2f\t%.2f\n", $sseArray[$i], $entry, $sse, $sum1sq;
    #}# Debug code end
  }
  #printf "bestSse = %.2f\t bestIndex=%d\n", $bestSse, $bestIndex;
  $m1 = &mean($data, $start, $bestIndex);
  $m2 = &mean($data, $bestIndex + 1, $end);
  my $thr = ($m1 + $m2)/2.0;
  #{# Debug code start
  #  $sum1sq = &mse($data, $start, $bestIndex);
  #  $sum2sq = &mse($data, $bestIndex + 1, $end);
  #  $sse = $sum1sq + $sum2sq;
  #  if ($sse != $bestSse) {
  #    print STDERR "SSE calculation is wrong :", $sse, ":", $bestSse, "\n";
  #  }
  #}# Debug code end
  my $label = 0;
  if ($m1 < $m2) {
    $label = 1;
  }
  else {
    $label = 2;
  }
  my $statistic = 0;
  if ($bestSse > 0) {
    if ($count > 4) {
      $statistic = ($sstot - $bestSse)/3/($bestSse/($count - 4));
    }
    else {
      $statistic = ($sstot - $bestSse)/2/$bestSse;
    }
  }
  return [$bestIndex, $bestSse, $sstot, $statistic, $m1, $m2, $thr, $label];
}

sub getStepMinerThr {
  my ($data, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$data) - 1;
  }
  my @array;
  for ($i = $start; $i <= $end; $i++) {
    next if (!defined $data->[$i] || $data->[$i] eq "");
    push @array, $data->[$i];
  }
  my @sorted = sort { $a <=> $b } @array;
  return &fitStep(\@sorted, 0, $#sorted);
}


sub printThr {
  my ($file, $start, $end, $gap) = @_;
  my $fh;
  open($fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  if ($end > $#headers) {
    $end = $#headers;
  }
  my $index = 0;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $res = &getStepMinerThr(\@list, $start, $end);
    print join("\t", $list[0], $res->[6], $res->[3], $res->[6]-$gap,
        $res->[6]+$gap), "\n";
    if ( ($index % 1000) == 0 ) {
      print STDERR $index, "\n";
    }
    $index++;
  }
  close($fh);
}

sub convertpcl {
  my ($file, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  my $fh;
  open($fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  if (!defined $end || $end > $#headers) {
    $end = $#headers;
  }
  print join("\t", @headers[0, 1]), "\tGWEIGHT\t", join("\t", @headers[$start .. $end]), "\n";
  print join("\t", "EWEIGHT", ""), "\t1\t", join("\t", map { 1 } @headers[$start .. $end]), "\n";
 
  my $index = 0;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
  print join("\t", @list[0, 1]), "\t1\t", join("\t", @list[$start .. $end]), "\n";
    if ( ($index % 1000) == 0 ) {
      print STDERR $index, "\n";
    }
    $index++;
  }
  close($fh);
}

sub printDynamicRangeInfo4 {
    my ($pclfile, $phfile, $phid) = @_;
    my $ph = &getPh($phfile, $phid);
    my $indexHash = {};
    for (my $i = 3; $i < scalar(@$ph); $i++) {
        if ($ph->[$i] == 1) {
            $indexHash->{$i} = 1;
        }
    }
    my $fh1;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @eweights = split(/\t/, $eweight);
    my @keys = sort keys(%{$indexHash});
    print join("\t", $headers[0], "0.5\%", "1\%", "2.5\%", "97.5\%", "99\%",
            "99.5\%", "95\% Dynamic-Range", "98\% Dynamic-Range",
            "99\% Dynamic-Range"), "\n";
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } @keys;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my @values = (0.005, 0.01, 0.025, 0.975, 0.99, 0.995);
      my @index = map { int($_*$#sortedExpr) } @values;
      my @data = map { $sortedExpr[$_] } @index;
      print join("\t", $list[0], @data,
          $data[3] - $data[2], $data[4] - $data[1], $data[5] - $data[0]), "\n";
    }
    close($fh1);
}

sub printDynamicRangeInfo5 {
    my ($pclfile, $thrfile, $thr2file, $phfile, $phid) = @_;
    my $maxPoints = 20; # -100 to 1, 0, 1, to 100
    my $gmin = 1.0;
    my $gmax = 16.0;
    my $gstep = 0.2;
    my $gnum = int(($gmax - $gmin)/$gstep);
    my $ph = &getPh($phfile, $phid);
    my $indexHash = {};
    for (my $i = 3; $i < scalar(@$ph); $i++) {
        if ($ph->[$i] == 1) {
            $indexHash->{$i} = 1;
        }
    }
    my $fh1;
    my $fh2;
    my $fh3;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    open($fh2, "<$thrfile") || die "Can't open $thrfile";
    open($fh3, "<$thr2file") || die "Can't open $thr2file";
    <$fh2>; <$fh3>;
    my $header=<$fh1>; 
    my $eweight=<$fh1>;
    $header =~ s/[\r\n]//g;
    $eweight =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    my @eweights = split(/\t/, $eweight);
    my @keys = sort keys(%{$indexHash});
    my @h1 = map { sprintf("Perc -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
    my @h2 = map { sprintf("Perc %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
    my @h3 = map { sprintf("Expr -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
    my @h4 = map { sprintf("Expr %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
    my @h5 = map { sprintf("Abs %.2f", $_*$gstep + $gmin) } 0 .. $gnum;
    print join("\t", $headers[0], "min");
    print "\t", join("\t", "0.5\%", "1\%", "2.5\%", "97.5\%", "99\%", "99.5\%");
    print "\t", join("\t", "max", "mean", "sd", "thr", "thrNum", "pval");
    print "\t", join("\t", @h1), "\t", join("\t", @h2);
    print "\t", join("\t", @h3), "\t", join("\t", @h4);
    print "\t", join("\t", @h5), "\n";
    my $lineno = 0;
    while (my $line = <$fh1>) {
      $lineno++;
      if (($lineno % 100) == 0) {
        print STDERR $lineno, "\r";
      }
      $line =~ s/[\r\n]//g;
      my $thrLine = <$fh2>;
      my ($n, $thr, $rest) = split(/\t/, $thrLine);    
      my $thr2Line = <$fh3>; $thr2Line =~ s/[\r\n]//g;
      my ($n, $pval, $rest) = split(/\t/, $thr2Line);    
      my @list = split(/\t/, $line);
      my @expr = map { $list[$_] } @keys;
      my @sortedExpr = sort {$a <=> $b} @expr;
      my ($min, $max, $mean, $sd, $num, $thrNum) = &getRangeInfo(\@expr, $thr);
      my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
      my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
      # Debug prints
      #print scalar(@sortedExprBelowThr), "\n";
      #print scalar(@sortedExprAboveThr), "\n";
      #print $thr, "\n";
      #print &BinarySearch([ 1, 2, 3, 4, 4, 5, 5, 6, 8, 9, 10], 5), "\n";
      print join("\t", $list[0], $min);
      my @values = (0.005, 0.01, 0.025, 0.975, 0.99, 0.995);
      my @indices = map { int($_*$#sortedExpr) } @values;
      my @data = map { $sortedExpr[$_] } @indices;
      print "\t", join("\t", @data);
      print "\t", join("\t", $max, $mean, $sd, $thr, $num - $thrNum, $pval);
      for (my $i = 0; $i < $maxPoints; $i++) {
        if (scalar(@sortedExprBelowThr) != 0) {
          my $index = int($i*1.0*$#sortedExprBelowThr/$maxPoints);
          my $expr = $sortedExprBelowThr[$index];
          print "\t", sprintf("%.2f", $expr);
        }
        else {
          print "\t", 0;
        }
      }
      for (my $i = 0; $i <= $maxPoints; $i++) {
        if (scalar(@sortedExprAboveThr) != 0) {
          my $index = int($i*1.0*$#sortedExprAboveThr/$maxPoints);
          my $expr = $sortedExprAboveThr[$index];
          if ($i == 0) {
            $expr = $thr;
          }
          print "\t", sprintf("%.2f", $expr);
        }
        else {
          print "\t", 0;
        }
      }
      my $min_below = $sortedExprBelowThr[0];
      my $max_below = $sortedExprBelowThr[$#sortedExprBelowThr];
      my $last_index = 0;
      for (my $i = 0; $i < $maxPoints; $i++) {
        if (scalar(@sortedExprBelowThr) != 0) {
          my $expr = $min_below + $i*($max_below - $min_below)/$maxPoints;
          my $index = &BinarySearch(\@sortedExprBelowThr, $expr);
          print "\t", sprintf("%d", $index - $last_index);
          $last_index = $index;
        }
        else {
          print "\t", 1;
        }
      }
      my $min_above = $sortedExprAboveThr[0];
      my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
      for (my $i = 0; $i <= $maxPoints; $i++) {
        if (scalar(@sortedExprAboveThr) != 0) {
          my $expr = $min_above + $i*($max_above - $min_above)/$maxPoints;
          my $index = &BinarySearch(\@sortedExprAboveThr, $expr);
          my $new_index = $index + scalar(@sortedExprBelowThr);
          print "\t", sprintf("%d", $new_index - $last_index);
          $last_index = $new_index;
        }
        else {
          print "\t", 1;
        }
      }
      my $last_index = 0;
      for (my $i = 0; $i <= $gnum; $i++) {
        if (scalar(@sortedExpr) != 0) {
          my $expr = $gmin + $i*$gstep;
          my $index = &BinarySearch(\@sortedExpr, $expr);
          if ($index < 0) {
            $index = 0;
          }
          print "\t", sprintf("%d", $index - $last_index);
          $last_index = $index;
        }
        else {
          print "\t", 0;
        }
      }
      print "\n";
    }
    close($fh1);
    close($fh2);
    close($fh3);
}

sub getMd5 {
  my ($filename) = @_;
  my $fh1;
  if ($filename =~ /gz$/) {
    open($fh1, "zcat '$filename'| md5sum |") || die "Can't open pipe $filename\n";
  }
  else {
    open($fh1, "md5sum '$filename'|") || die "Can't open pipe $filename\n";
  }
  my $line = <$fh1>;
  my ($num, $f) = split(/\s+/, $line);
  close($fh1);
  return $num;
}

sub getMd5Sum {
  my ($file, $igfile) = @_;
  my $ighash = {};
  if (defined $igfile && -e $igfile) {
    my $fh;
    open($fh, "<$igfile") || die "Can't open $igfile\n";
    while (<$fh>) {
      s/[\r\n]//g;
      my $id = &getMd5($_);
      $ighash->{$id} = $_;
    }
    close($fh);
  }

  my $hash1 = {};
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $index = 0;
  while (<$fh>) {
    s/[\r\n]//g;
    my $filename = $_;
    my $id = &getMd5($filename);
    if (!defined $ighash->{$id}) {
      if (!defined $hash1->{$id}) {
        $hash1->{$id} = $filename;
        print "$filename\t$id\n";
      }
      else {
        print "$filename\t$id\tdup\t". $hash1->{$id} ."\n";
      }
    }
    if ( ($index % 1000) == 0) {
      print STDERR "$index\n";
    }
    $index++;
  }
  close($fh);
}

sub mergeMd5Sum {
  my @files = @_;
  my $hash1 = {};
  my $index = 0;
  foreach my $file (@files) {
    open(my $fh, "<$file") || die "Can't open $file\n";
    while (<$fh>) {
      s/[\r\n]//g;
      my @list = split("\t");
      my $filename = $list[0];
      my $id = $list[1];
      if ($list[2] ne "dup") {
        if (!defined $hash1->{$id}) {
          $hash1->{$id} = $filename;
          print "$filename\t$id\n";
        }
        else {
          print "$filename\t$id\tdup\t". $hash1->{$id} ."\n";
        }
      }
      else {
        print join("\t", @list), "\n";
      }
      if ( ($index % 1000) == 0) {
        print STDERR "$index\n";
      }
      $index++;
    }
    close($fh);
  }
}

sub selectNoDup {
  my ($file, $igfile, $rest) = @_;
  my $hash = {};
  if (defined $igfile) {
    open(my $fh, "<$igfile") || die "Can't open $igfile\n";
    while (<$fh>) {
      s/[\r\n]//g;
      my @list = split("\t");
      if ($list[2] ne "dup") {
        $hash->{$list[1]} = $list[0];
      }
    }
    close($fh);
  }
  open(my $fh, "<$file") || die "Can't open $file\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    next if (defined $hash->{$list[1]});
    if ($list[2] ne "dup") {
        print $list[0], "\n";
    }
  }
  close($fh);
}

sub selectPlatform {
  my ($file, $type, $listfile, $igfile) = @_;
  my $ighash = {};
  open(my $fh, "<$igfile");
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    $ighash->{$list[0]} = 1;
  }
  close($fh);
  my $hash = {};
  open(my $fh, "<$listfile");
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    $hash->{$list[0]} = 1;
  }
  close($fh);
  open(my $fh, "<$file") || die "Can't open $file\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    next if (defined $listfile && !defined $hash->{$list[1]});
    next if (defined $igfile && defined $ighash->{$list[1]});
    if ($list[0] eq $type) {
        print $list[1], "\n";
    }
  }
  close($fh);
}

sub getWString {
  my ($fh, $rest) = @_;
  my $buffer;
  read($fh, $buffer, 4); my $n = unpack 'N', $buffer;
  read($fh, $buffer, $n*2); 
  my @arr = split(//, $buffer);
  $buffer = join("", map {$arr[$_]} map {$_*2+1} 0 .. $n);
  my $str = unpack 'a*', $buffer;
  return ($n, $str);
}

sub getString {
  my ($fh, $rest) = @_;
  my $buffer;
  read($fh, $buffer, 4); my $n = unpack 'N', $buffer;
  read($fh, $buffer, $n); my $str = unpack 'a*', $buffer;
  return ($n, $str);
}

sub getCelInfo {
  my $file = shift;
  my $buffer;
  my $res;
  my $fh;
  if ($file =~ /.gz$/) {
    open($fh, "zcat '$file' |") || die "Can't open zcat $file\n";
  }
  else {
    open($fh, "<$file") || die "Can't open $file\n";
  }
  binmode($fh);
  read($fh, $buffer, 1);
  my $magic = unpack 'C', $buffer;
  if ($magic == 59) {
    read($fh, $buffer, 1); my $version = unpack 'c', $buffer;
    read($fh, $buffer, 4); my $numDataGroup = unpack 'N', $buffer;
    read($fh, $buffer, 4); my $posDataGroup = unpack 'N', $buffer;
    #print STDERR "AGCC $version - $numDataGroup - $posDataGroup\n";
    my ($n, $str) = &getString($fh);
    #print STDERR "$n: Id $str :-\n";
    my ($n, $str) = &getString($fh);
    #print STDERR "$n: Guid $str :-\n";
    my ($n, $str) = &getWString($fh);
    #print STDERR "$n: time $str :-\n";
    my ($n, $str) = &getWString($fh);
    #print STDERR "$n: locale $str :-\n";
    read($fh, $buffer, 4); my $n = unpack 'N', $buffer;
    #print STDERR "Num par $n :-\n";
    for (my $i = 0; $i < $n; $i++) {
      last if ($i > 100);
      my ($name, $value, $type);
      my ($n, $name) = &getWString($fh);
      #print STDERR "$n: Name $name :-", length($name), "\n";
      my ($n, $value) = &getString($fh);
      $value =~ s/\000//g;
      #print STDERR "$n: Value $value :-\n";
      my ($n, $type) = &getWString($fh);
      #print STDERR "$n: Type $type :-\n";
      if ($name =~ /affymetrix-partial-dat-header/ ||
          $name =~ /affymetrix-array-type/) {
        return $value;
      }
    }
  }
  while (<$fh>) {
    if (/DatHeader=/) {
      if (/([^\s]*).1sq/) {
        $res = $1;
      }
    }
  }
  close($fh);
  return $res;
}

sub printCELinfo {
  my ($file, $info) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $f = $list[0];
    my $info = &getCelInfo($f);
    print join("\t", $info, $f), "\n";
  }
  close($fh);
}

sub getCelVersion {
  my $file = shift;
  my $buffer;
  my $res;
  my $fh;
  if ($file =~ /.gz$/) {
    open($fh, "zcat '$file' |") || die "Can't open zcat $file\n";
  }
  else {
    open($fh, "<$file") || die "Can't open $file\n";
  }
  binmode($fh);
  read($fh, $buffer, 1);
  my $magic = unpack 'C', $buffer;
  if ($magic == 59) {
    return "AGCC 1";
  }
  elsif ($magic == 64) {
    return 4;
  }
  elsif ($magic == 91) {
    return 3;
  }
  close($fh);
  return -1;
}

sub printCELversion {
  my ($file, $info) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $f = $list[0];
    my $info = &getCelVersion($f);
    print join("\t", $info, $f), "\n";
  }
  close($fh);
}

sub selectNew {
  my ($file, @files) = @_;
  my $hash = {};
  my $index = 0;
  foreach my $f (@files) {
    open(my $fh, "<$f") || die "Can't open $f\n";
    while (<$fh>) {
      s/[\r\n]//g;
      my @list = split("\t");
      my $filename = $list[0];
      $hash->{$filename} = 1;
      if ( ($index % 1000) == 0) {
        print STDERR "$index\n";
      }
      $index++;
    }
    close($fh);
  }
  open(my $fh, "<$file") || die "Can't open $file\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    next if (defined $hash->{$list[0]});
    print $list[0], "\n";
    $hash->{$list[0]} = 1;
  }
  close($fh);
}


