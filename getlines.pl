#! /usr/bin/perl

# script created by Sorel Fitz-Gibbon
# to use the script, 

if ($ARGV[1] eq undef) {
   die "Give <list of names> and <file of lines with name> . MUST EDIT TO SPECIFY WHERE TO FIND NAME. Output will be $ARGV[0].lines \n";
}

open(IN,$ARGV[0]) || die;
while(<IN>) {
#   /^(\S+)/;
   chomp;
   $keepers{$_} = 1;
}

open(IN,"$ARGV[1]") || die;
open(OUT,">$ARGV[0].lines");
while(<IN>) {
      chomp;
   unless (/^chr/) {
      /^(m01oak[^-]+)-/ || die;
      $check = $1;
     if ($keepers{$check}) {
         delete $keepers{$check};
         print OUT "$_\n";
     }
   }
}

@leftovers = (keys %keepers);
if ($#leftovers == -1) {
   print "all input patterns found\n";
} else {
   $notfound = $#leftovers+1;
   print "$notfound input patterns not found\n";
}
unless ($#leftovers < 0 ) {
   print "These lines not found...\n";
   for $line (@leftovers) {
      print "$line\n";
   }
}
