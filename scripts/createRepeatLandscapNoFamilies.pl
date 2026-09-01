#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) createRepeatLandscape.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Juan Caballero     jcaballero@systemsbiology.org
##  Description:
##      Based on earlier work by Juan Caballero
##
#******************************************************************************
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************************************************
#
# ChangeLog
#
#     $Log$
#
###############################################################################
#
# To Do:
#

=head1 NAME

createRepeatLandscape.pl - Create a Repeat Landscape graph

=head1 SYNOPSIS

  createRepeatLandscape.pl [-version] -div *.divsum [-t "graph title"]
                           [-j] -g # | -twoBit <twoBitFile>

=head1 DESCRIPTION

  Create a Repeat Landscape graph using the divergence summary data
  generated with the calcDivergenceFromAlign.pl script.


=head1 EXAMPLES

  Older ( pre 4.0.4 ) RepeatMasker dataset:

     ./calcDivergenceFromAlign.pl -s example.divsum -a example_with_div.align 
                                  example.align.gz
 
     This creates an additional file "example_with_div.align" which contains
     the added Kimura divergence field after each alignment.

     ./createRepeatLandscape.pl -div example.divsum > 
                                /home/user/public_html/example.html 


  On newer RepeatMasker dataset that already contains the Kimura divergence
  line following each alignment:

     ./calcDivergenceFromAlign.pl -s example.divsum example.align.gz

     ./createRepeatLandscape.pl -div example.divsum > 
                                /home/user/public_html/example.html 


The options are:

=over 4

=item -version

Displays the version of the program

=item -div <file>

The divergence summary file created with the calcDivergenceFromAlign.pl 
script.

=item -g #

Set the genome size used in percentage calculations.

=item -twoBit <filename>

Get the genome size directly from the sequence file 
( excluding Ns ).  This option requires that the
UCSC utility "twoBitInfo" is in your path.

=item -j 

Output javascript only and not a fully constructed HTML page.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2013-2014 Robert Hubley, Institute for Systems Biology

=head1 AUTHORS

Robert Hubley <rhubley@systemsbiology.org>

Juan Caballero <jcaballero@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Data::Dumper;
use Getopt::Long;

#
# Version
#
#  This is a neat trick.  CVS allows you to tag
#  files in a repository ( i.e. cvs tag "2003/12/03" ).
#  If you check out that release into a new
#  directory with "cvs co -r "2003/12/03" it will
#  place this string into the $Name:  $ space below
#  automatically.  This will help us keep track
#  of which release we are using.  If we simply
#  check out the code as "cvs co Program" the
#  $Name:  $ macro will be blank so we should default
#  to what the ID tag for this file contains.
#
my $CVSNameTag = '$Name:  $';
my $CVSIdTag   =
    '$Id: createRepeatLandscape.pl,v 1.14 2017/02/01 21:01:57 rhubley Exp $';
my $Version = $CVSNameTag;
$Version = $CVSIdTag if ( $Version eq "" );

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',      # print out the version and exit
                    '-div=s',
                    '-maxScale=s',
                    '-minCount=s',
                    '-g=s',
                    '-twoBit=s',
                    '-t=s',
                    '-j'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

if ( !-s $options{'div'} ) {
  print "\n\nMissing '-div' option or the file is missing\n";
  usage();
}

my $kout  = $options{'div'};
my $gsize = 0;
if ( $options{'g'} ) {
  $gsize = $options{'g'};
}
elsif ( $options{'twoBit'} ) {
  open IN, "twoBitInfo -noNs $options{'twoBit'} stdout|"
      or die
      "Could not open twoBitInfo: twoBitInfo -noNs $options{'twoBit'} stdout\n";
  while ( <IN> ) {
    if ( /^\S+\s+(\d+)/ ) {
      $gsize += $1;
    }
  }
  close IN;
}
else {
  print "\n\nMust provide either the -g or the -twoBit flags!\n\n";
  usage();
}
my $out     = undef;
my $help    = undef;
my $verbose = undef;

# Main variables
my %size;
my %divData;
my %dataClasses;
my %div;
my %sumDiv;

# order is IMPORTANT for plotting

my $graphLabels = [
                    [ 'Unknown',        '#999999' ],
                    [ 'Other',          '#4D4D4D' ],
                    [ 'DNA/Academ',     '#FF0000' ],
                    [ 'DNA/CMC',        '#FF200B' ],
                    [ 'DNA/Crypton',    '#FF3115' ],
                    [ 'DNA/Ginger',     '#FF3D1E' ],
                    [ 'DNA/Harbinger',  '#FF4825' ],
                    [ 'DNA/hAT',        '#FF512D' ],
                    [ 'DNA/Kolobok',    '#FF5A34' ],
                    [ 'DNA/Maverick',   '#FF623B' ],
                    [ 'DNA',            '#FF6A42' ],
                    [ 'DNA/Merlin',     '#FF7149' ],
                    [ 'DNA/MULE',       '#FF7850' ],
                    [ 'DNA/P',          '#FF7F57' ],
                    [ 'DNA/PiggyBac',   '#FF865E' ],
                    [ 'DNA/Sola',       '#FF8D65' ],
                    [ 'DNA/TcMar',      '#FF936C' ],
                    [ 'DNA/Transib',    '#FF9972' ],
                    [ 'DNA/Zator',      '#FF9F79' ],
                    [ 'DNA/Dada',       '#FFCFBC' ],
                    [ 'DIRS',    '#FF00FF' ],
                    [ 'LTR/DIRS',       '#006400' ],
                    [ 'LTR/Ngaro',      '#197214' ],
                    [ 'LTR/Pao',        '#2A8024' ],
                    [ 'LTR/Copia',      '#3A8F33' ],
                    [ 'LTR/Gypsy',      '#489E42' ],
                    [ 'LTR/ERVL',       '#57AE51' ],
                    [ 'LTR',            '#65BD61' ],
                    [ 'LTR/ERV1',       '#73CD70' ],
                    [ 'LTR/ERV',        '#81DD80' ],
                    [ 'LTR/ERVK',       '#90ED90' ],
                    [ 'LINE/L1',        '#00008B' ],
                    [ 'LINE',           '#251792' ],
                    [ 'LINE/RTE',       '#38299A' ],
                    [ 'LINE/CR1',       '#483AA2' ],
                    [ 'LINE/Rex-Babar', '#554BAA' ],
                    [ 'LINE/L2',        '#625CB1' ],
                    [ 'LINE/Proto2',    '#6E6DB9' ],
                    [ 'LINE/LOA',       '#797EC0' ],
                    [ 'LINE/R1',        '#848FC8' ],
                    [ 'LINE/Jockey-I',  '#8FA1CF' ],
                    [ 'LINE/Dong-R4',   '#99B3D7' ],
                    [ 'LINE/R2',        '#A3C5DE' ],
                    [ 'LINE/CRE',       '#C1D9FF' ],
                    [ 'PLE',            '#ACD8E5' ],
                    [ 'Retroposon/SVA', '#FF4D4D' ],
                    [ 'SINE',           '#9F1FF0' ],
                    [ 'SINE/5S',        '#A637F1' ],
                    [ 'SINE/7SL',       '#AD49F2' ],
                    [ 'SINE/Alu',       '#B358F3' ],
                    [ 'SINE/tRNA',      '#B966F4' ],
                    [ 'SINE/tRNA-Alu',  '#BF74F4' ],
                    [ 'SINE/tRNA-RTE',  '#C481F5' ],
                    [ 'SINE/RTE',       '#C98EF6' ],
                    [ 'SINE/Deu',       '#CE9BF7' ],
                    [ 'SINE/tRNA-V',    '#D3A7F7' ],
                    [ 'SINE/MIR',       '#D7B4F8' ],
                    [ 'SINE/U',         '#DFCDF9' ],
                    [ 'SINE/tRNA-7SL',  '#E2D9F9' ],
                    [ 'SINE/tRNA-CR1',  '#E5E5F9' ]
];

my $extraChartLabels = [
                         [ 'Simple_repeat',  '#FFFF66' ],
                         [ 'Satellite',      '#CCCC52' ],
                         [ 'Structural_RNA', '#99993D' ],
                         [ 'Low_complexity', '#666629' ]
];

my $rec    = 0;
my $masked = 0;
my @rep    = ();
print STDERR "Parsing $kout\n";
open F, "$kout" or die "cannot read $kout\n";
my $inCoverageSection = 0;
my $inClassSection    = 0;
my %classTotalBP      = ();
while ( <F> ) {
  chomp;
  my @fields = split( /\s+/, $_ );

  if ( /^Coverage for each repeat/ ) {
    $inClassSection    = 0;
    $inCoverageSection = 1;
    next;
  }

  #if ( $gsize == 0 && /^Genome Size\s+=\s+(\d+)/ ) {
  #  $gsize = $1;
  #}
  if ( /^Class\tRepeat/ ) {
    $inClassSection = 1;
    next;
  }

  if ( $inClassSection && /^\S+\t\S+/ ) {
    next if ( /^[-\s]+$/ );
    my @fields = split( /\t/ );
    if ( @fields == 5 ) {

      # Class	Repeat	absLen	wellCharLen	Kimura%
      my $class = fixName( $fields[ 0 ] );
      $classTotalBP{$class} += $fields[ 2 ];
    }
    else {
      die
"\nThis divergence summary file ( $options{'div'} ) was created with an\nolder version "
          . "of calcDivergenceFromAlign.pl.  Please regenerate this file\nusing the current version of "
          . "calcDivergenceFromAlign.pl and then rerun this script.\n\n";
    }
  }
  elsif ( $inCoverageSection && /^Div\s+\S/ ) {

    # Get rid of "Div"
    shift @fields;
    foreach my $class ( @fields ) {
      $class = fixName( $class );
      $dataClasses{$class} = 1;
      push @rep, $class;
    }
    $rec = 1;
  }
  elsif ( $rec == 1 ) {
    my $div = shift @fields;
    $div{$div} = 1;
    for ( my $i = 0 ; $i <= $#fields ; $i++ ) {
      my $per = 100 * $fields[ $i ] / $gsize;
      my $r   = $rep[ $i ];
      $divData{$r}{$div} += $per;
      $sumDiv{$r}        += $per;
    }
  }
}
close F;

my $minCount = 0;
if ( $options{'minCount'} ){
  $minCount = $options{'minCount'};
}

## Double check for new repeat classes we don't have in our
## fixed graph label list.
# Create a lookup hash of the graph labels
my %graphLabelHash = ();
foreach my $label ( @{$graphLabels} ) {
  $graphLabelHash{ $label->[ 0 ] }++;
}
foreach my $key ( keys %dataClasses ) {
  if ( !defined $graphLabelHash{$key}
       && $key !~ /^(Satellite|Segmental|Structural_RNA|Simple_repeat)$/ )
  {
    warn "$key will not be graphed!\n";
  }
}

my $colors               = "";
my $pieColors            = "slices: {";
my $idx                  = 0;
my @landscapeGraphLabels = ();
my @div                  = sort { $a <=> $b } keys %div;
foreach my $label ( @{$graphLabels} ) {
  if ( defined $sumDiv{ $label->[ 0 ] }
       && $sumDiv{ $label->[ 0 ] } > 0 
       && defined $classTotalBP{ $label->[ 0 ] }
       && $classTotalBP{ $label->[ 0 ] } > $minCount )
  {
    push @landscapeGraphLabels, $label->[ 0 ];
    $colors .= "\"$label->[1]\", ";
  }
}

foreach my $label ( reverse( @{$graphLabels}, @{$extraChartLabels} ) ) {
  if ( defined $classTotalBP{ $label->[ 0 ] } && $classTotalBP{ $label->[ 0 ] } > $minCount ) {
    $pieColors .= "$idx: { color: \"$label->[1]\" }, ";
    $idx++;
  }
}

# For unmasked part of pie chart
$pieColors .= "$idx: { color: \"black\" }} ";
$colors =~ s/, $//;

#
# div  class1 class2..
# frac
#
my %pdata = ();
foreach my $div ( @div ) {
  next if $div > 50;
  my @tmp = ();
  foreach my $rep ( @landscapeGraphLabels ) {
    my $per = 0;
    $per = $divData{$rep}{$div} if ( defined $divData{$rep}{$div} );
    push @tmp, $per;
  }
  $pdata{$div} = join( ', ', @tmp );
}

# Default parameters
my $title = 'Interspersed Repeat Landscape';
$title = $options{'t'} if ( $options{'t'} );

my @genomeFractionChartLabels = @landscapeGraphLabels;
foreach my $label ( @{$extraChartLabels} ) {
  if ( defined $classTotalBP{ $label->[ 0 ] }
       && $classTotalBP{ $label->[ 0 ] } > $minCount )
  {
    push @genomeFractionChartLabels, $label->[ 0 ];
  }
}

# Make it safe for the web
$title =~ s/_/ /g;

unless ( $options{'j'} ) {
  print <<_HTML_HEADER_
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
  <head>
    <title>Interspersed Repeat Landscape</title>
    <script type="text/javascript" src="https://www.google.com/jsapi"></script>
    <script type="text/javascript">
_HTML_HEADER_
      ;
}

print <<_HEADER_
      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = new google.visualization.DataTable();
        data.addColumn('string', 'Divergence');
_HEADER_
    ;

foreach my $name ( @landscapeGraphLabels ) {
  print "          data.addColumn('number', '$name');\n";
}

print "          data.addRows(\[\n";
foreach my $div ( sort { $a <=> $b } keys %pdata ) {
  my $dd = $pdata{$div};
  print "          \[\'$div\', $dd\],\n";
}
print "        ]);\n";

print "          var pieData = new google.visualization.DataTable();\n";
print "          pieData.addColumn('string', 'class');\n";
print "          pieData.addColumn('number', 'bp');\n";
print "          pieData.addRows(\[\n";
foreach my $rClass ( reverse @genomeFractionChartLabels ) {
  print "         \[\'$rClass\', $classTotalBP{$rClass}\],\n";
  $masked += $classTotalBP{$rClass};
}
my $unmasked = $gsize - $masked;
print "              \[\'Unmasked\',$unmasked\],\n";
print "        ]);\n";

print <<_PIE_
        var pieOptions = {
          legend: 'none',
          title: 'Genome Fraction',
          $pieColors
        };
        var pieChart = new google.visualization.PieChart(document.getElementById('pie_chart_div'));
        pieChart.draw(pieData, pieOptions);
_PIE_
    ;

print <<_TAIL_
        var options = {
          animation: {duration: 10},
          title: '$title',
          hAxis: {title: 'Kimura substitution level (CpG adjusted)', showTextEvery: 5},
          vAxis: {title: 'percent of genome'},
          isStacked: 1,
          colors: [$colors]
        };
        var chart = new google.visualization.ColumnChart(document.getElementById('chart_div'));
        chart.draw(data, options);
      }
_TAIL_
    ;

unless ( $options{'j'} ) {
  print <<_HTML_TAIL_
      </script>
    </head>
    <body bgcolor="#ffffff" text="#000000" link="#525d76">
    <font size="+3">Interspersed Repeat Landscape</font>
   <hr noshade="noshade" size="1">
   <div id="outer" style="width:9999px">
   <div id="chart_div" style="width: 1000px; height: 600px;float:left"></div>
   <div id="pie_chart_div" style="width: 500px; height: 600px;float:left"></div>
   </div>
   <p>&#169; RepeatMasker.org</p>  
   </body>
</html>
_HTML_TAIL_
      ;
}

sub fixName {
  my $r = shift;
  $r =~ s/\?//g;

# Complete simplified TE class mapping (major groups only)

my %nameMap = (
    # DNA transposons (all → "DNA")
    "DNA/Chompy"           => "DNA",
    "DNA/CMC-Chapaev"      => "DNA",
    "DNA/CMC-Chapaev-3"    => "DNA",
    "DNA/CMC-EnSpm"        => "DNA",
    "DNA/CMC-Transib"      => "DNA",
    "DNA/En-Spm"           => "DNA",
    "DNA/Crypton-3"        => "DNA",
    "DNA/PIF-Harbinger"    => "DNA",
    "DNA/PIF-ISL2EU"       => "DNA",
    "DNA/Tourist"          => "DNA",
    "DNA_PIF-Harbinger"    => "DNA",
    "PIF-Harbinger"        => "DNA",
    "DNA/AcHobo"           => "DNA",
    "DNA/Charlie"          => "DNA",
    "DNA/Chompy1"          => "DNA",
    "DNA/MER1_type"        => "DNA",
    "DNA/Tip100"           => "DNA",
    "DNA/hAT-Ac"           => "DNA",
    "DNA/hAT-Blackjack"    => "DNA",
    "DNA/hAT-Charlie"      => "DNA",
    "DNA/hAT-Tag1"         => "DNA",
    "DNA/hAT-Tip100"       => "DNA",
    "DNA/hAT-hATw"         => "DNA",
    "DNA/hAT-hobo"         => "DNA",
    "DNA/hAT_Tol2"         => "DNA",
    "DNA_hAT-Ac"           => "DNA",
    "DNA/hAT-Pegasus"      => "DNA",
    "DNA/Kolobok-IS4EU"    => "DNA",
    "DNA/Kolobok-T2"       => "DNA",
    "DNA/T2"               => "DNA",
    "DNA/MULE-MuDR"        => "DNA",
    "DNA/MULE-NOF"         => "DNA",
    "DNA/MuDR"             => "DNA",
    "DNA/piggyBac"         => "DNA",
    "DNA_PiggyBac"         => "DNA",
    "DNA_Piggybac"         => "DNA",
    "DNA/MER2_type"        => "DNA",
    "DNA/Mariner"          => "DNA",
    "DNA/Pogo"             => "DNA",
    "DNA/Stowaway"         => "DNA",
    "DNA/Tc1"              => "DNA",
    "DNA/Tc2"              => "DNA",
    "DNA/Tc4"              => "DNA",
    "DNA/TcMar-Fot1"       => "DNA",
    "DNA/TcMar-ISRm11"     => "DNA",
    "DNA/TcMar-Mariner"    => "DNA",
    "DNA/TcMar-Pogo"       => "DNA",
    "DNA/TcMar-Tc1"        => "DNA",
    "DNA/TcMar-Tc2"        => "DNA",
    "DNA/TcMar-Tigger"     => "DNA",
    "DNA/Tigger"           => "DNA",
    "DNA_Tc1-Mariner"      => "DNA",
    "DNA_Tc1Mariner"       => "DNA",
    "DNA_TcMar-Tc1"        => "DNA",
    "DNA/Zisupton"         => "DNA",
    "DNA/Helitron"         => "DNA",
    "DNA_Zator"            => "DNA",
    "DNA/CMC-3"            => "DNA",
    "DNA/CMC-Mirage"       => "DNA",
    "DNA/Dada"             => "DNA",
    "DNA/Kolobok-Hydra"    => "DNA",
    "DNA/MULE-F"           => "DNA",
    "DNA/TcMar-Stowaway"   => "DNA",
    "DNA/TcMar-Tc4"        => "DNA",
    "DNA/TcMar-m44"        => "DNA",
    "DNA/hAT-Pegasus"      => "DNA",
    "DNA/hAT-Tol2"         => "DNA",
    "DNA/hAT-hAT1"         => "DNA",
    "DNA/hAT-hAT5"         => "DNA",
    "DNA/hAT-hAT6"         => "DNA",
    "DNA/hAT-hAT19"        => "DNA",
    "DNA_CACTA"            => "DNA",
    "DNA_hAT"              => "DNA",
    "DNA_TcMar"            => "DNA",
    "DNA/PiggyBac"         => "DNA",
    "DNA/hAT"              => "DNA",
    

    # LTR retrotransposons (all → "LTR")
    "LTR/DIRS1"            => "DIRS",
    "LTR/DIRS"            => "DIRS",
    "LTR/ERV-Foamy"        => "LTR",
    "LTR/ERV-Lenti"        => "LTR",
    "LTR/ERVL-MaLR"        => "LTR",
    "LTR/Gypsy-Troyka"     => "LTR",
    "LTR/MaLR"             => "LTR",
    "LTR_Gypsy"            => "LTR",
    "LTR/ERV4"             => "LTR",
    "LTR/Copia(Xen1)"      => "LTR",
    "LTR/Ginger"           => "LTR",
    "DIRS_DIRS1"           => "DIRS",
    "LTR/ERV"              => "LTR",
    "LTR/ERV1"             => "LTR",
    "LTR/ERVL"             => "LTR",
    "LTR/Gypsy"            => "LTR",
    "LTR/Copia"            => "LTR",
    "LTR/Pao"              => "LTR",
    "LTR/Ngaro"            => "LTR",

    # LINE elements (all → "LINE")
    "LINE/CR1-Zenon"       => "LINE",
    "LINE/I"               => "LINE",
    "LINE/Jockey"          => "LINE",
    "LINE/L1-Tx1"          => "LINE",
    "LINE/R2-Hero"         => "LINE",
    "LINE/RTE-BovB"        => "LINE",
    "LINE/RTE-RTE"         => "LINE",
    "LINE/RTE-X"           => "LINE",
    "LINE/telomeric"       => "LINE",
    "LINE/Penlope"         => "LINE",
    "LINE/Penelope"        => "LINE",
    "LINE_CR1"             => "LINE",
    "LINE_L1"              => "LINE",
    "LINE_L2"              => "LINE",
    "LINE/CRE"             => "LINE",
    "LINE/Jockey-I-I"      => "LINE",
    "LINE/I-Jockey"        => "LINE",
    "LINE/R2-NeSL"         => "LINE",
    "LINE/L2"              => "LINE",
    "LINE/Rex-Babar"       => "LINE",
    "LINE/CR1"             => "LINE",
    "LINE/L1"              => "LINE",
    
    # SINE elements (all → "SINE")
    "SINE/B2"              => "SINE",
    "SINE/B4"              => "SINE",
    "SINE/BovA"            => "SINE",
    "SINE/C"               => "SINE",
    "SINE/Core"            => "SINE",
    "SINE/ID"              => "SINE",
    "SINE/Lys"             => "SINE",
    "SINE/MERMAID"         => "SINE",
    "SINE/RTE-BovB"        => "SINE",
    "SINE/tRNA-Glu"        => "SINE",
    "SINE/tRNA-Lys"        => "SINE",
    "SINE/V"               => "SINE",
    "SINE/5S-Deu-L2"       => "SINE",
    "SINE/5S-Sauria-RTE"   => "SINE",
    "SINE/L2"              => "SINE",
    "SINE/Mermaid"         => "SINE",
    "SINE/R2"              => "SINE",
    "SINE/U"               => "SINE",
    "SINE/Core-RTE"        => "SINE",
    "SINE/tRNA-C"          => "SINE",
    "SINE/tRNA-Core-RTE"   => "SINE",
    "SINE/tRNA-Core"       => "SINE",
    "SINE/tRNA-Deu-CR1"    => "SINE",
    "SINE/tRNA-Deu-L2"     => "SINE",
    "SINE/tRNA-Deu"        => "SINE",
    "SINE/tRNA-Jockey"     => "SINE",
    "SINE/tRNA-L2"         => "SINE",
    "SINE/tRNA-Rex"        => "SINE",
    "SINE/tRNA-Sauria-L2"  => "SINE",
    "SINE/tRNA-Sauria-RTE" => "SINE",
    "SINE/tRNA-Sauria"     => "SINE",
    "SINE/Sauria"          => "SINE",
    "SINE/tRNA-V-Core-L2"  => "SINE",
    "SINE/tRNAore-RTE"     => "SINE",
    "SINE/tRNAore"         => "SINE",

    # Other categories
    "RC/Helitron"          => "RC",
    "RC/Helitron-2"        => "RC",
    "PLE/Hydra"            => "PLE",
    "PLE/Naiad"            => "PLE",
    "PLE/Athena"           => "PLE",
    "PLE/Poseidon"         => "PLE",
    "PLE/Nematis"          => "PLE",
    "PLE/Neptune"          => "PLE",
    "PLE/Coprina"          => "PLE",
    "Unknown/Y-chromosome" => "Unknown",
    "Retroposon"           => "Other",
    "Other/Composite"      => "Other",
    "Satellite/centromeric"=> "Satellite",
    "Satellite/telomeric"  => "Satellite",
    "Satellite/acromeric"  => "Satellite",
    "tRNA"                 => "Structural_RNA",
    "scRNA"                => "Structural_RNA",
    "RNA"                  => "Structural_RNA",
    "snRNA"                => "Structural_RNA",
    "rRNA"                 => "Structural_RNA"
);



  if ( $nameMap{$r} ) {
    $r = $nameMap{$r};
  }
  return $r;
}

1;
