#!/usr/bin/perl -w

# Extract the sun or spacecraft position from a list of
# cube files. The list of cube files must be one file per line.
# The output goes to a text file.

# This script assumes that the environmental variables
# ISISROOT and ISIS3DATA have been set up, and the
# ISIS executable campt is in the path.

use strict;

if (scalar(@ARGV) < 3){
  print_usage();
  exit(1);
}

my $type = $ARGV[0];
if ($type ne "sun" && $type ne "spacecraft"){
  print_usage();
  exit(1);
}
  
my $cubesList = $ARGV[1];
open(FILE, "<$cubesList") || die "ERROR: Cannot find file: $cubesList\n";
my @cubes = <FILE>;
close(FILE);

my $outputPositionList = $ARGV[2];
open(FILE, ">$outputPositionList") || die "ERROR: Cannot write to file $outputPositionList\n";

if ( !exists $ENV{'ISISROOT'} || !exists $ENV{'ISIS3DATA'} ){
  print "ERROR: The environmental variables ISISROOT and ISIS3DATA need to be set.\n";
  exit(1);
}

my $campt_exe = qx(which campt 2>&1); $campt_exe =~ s/\n//g;
if ($campt_exe =~ /(not found|no campt)/){
  print "ERROR: The ISIS executable campt was not found.\n";
  exit(1);
}

foreach my $file (@cubes){
  $file =~ s/\n//g;
  #print $file . "\n";
  my $cmd = "\$ISISROOT/scripts/isis3Startup.sh; $campt_exe from=" . $file;
  print $cmd . "\n";
  my $info = qx($cmd);
  my $pos    = extract_pos($type, $info);
  #print $pos . "\n";
  print FILE $pos . "\n";
}
close(FILE);

sub extract_pos{

  my $type    = shift;
  my $infoTxt = shift;
  my @info    = split("\n", $infoTxt);
  my $pos     = "";
  
  my $forcedCopy         = 0;
  my $forcedFilenameCopy = 0;
  foreach my $line (@info) {
    
    my @vals = split(' ', $line);
     
    if ($forcedFilenameCopy == 1) {
      $pos .= $vals[0] . " ";
      $forcedFilenameCopy = 0;
    }

    if (scalar(@vals) >= 3 && $vals[0] eq 'Filename') {
      $pos .= $vals[2];
      $forcedFilenameCopy = 1; 
    }

    if ($forcedCopy == 1) {
      $pos .= $vals[0] . "\n";
      $forcedCopy = 0;
    }
    if (scalar(@vals) >= 4 &&
        (
         ($type eq 'spacecraft' && $vals[0] eq 'SpacecraftPosition') ||
         ($type eq 'sun'        && $vals[0] eq 'SunPosition')
        )
        
       ) {
      $pos .= $vals[2] . " " . $vals[3] . " ";
      $forcedCopy = 1;
    }
         
  }

  # Convert to the format: ASXX-M-XXX x y z
  $pos =~ s/[,\(\)\s]+/ /g;
  $pos =~ s/^.*?(AS\d+-M-\d+).*? /$1 /g;
  $pos =~ s/\s*$//g;

  return $pos;
}

sub print_usage{
  print "Usage: $0 sun-or-spacecraft inputCubesList.txt outputPositionList.txt\n";
}
