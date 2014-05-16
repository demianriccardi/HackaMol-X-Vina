# Demian Riccardi May 13, 2014
#
# This script is for running virtual screens using Autodock Vina.  This is a work in
# progress!  Eventually, the functionality will be encapsulated in a class for screening.
# Obviously, this needs to be reworked to use a proper database!
#
# INPUT:
#   this script reads in a YAML configuration file passed on commandline:
#   $yaml->{in_json} is the path to the json file containing independent data
#   for each ligand. The first run, there will be no docking information present.
#   The json file will contain the ligand information (TMass, formula, BEST => {BE = 0},
#   etc. ). This script iterates (see JSON::XS on metacpan) through this json loaded by 
#   ligand and running the centers and receptors from the YAML configuration file.  
#
#   The set of centers are assumed to correspond to the set of receptors! 
#   e.g. two receptors with very different coordinates (translation) should have two
#        different sets of centers.
#
#   The num_modes will always be 1 for screens; i don't see the point (yet) of additional 
#   configs for virtual screens.  Need a different script to run with num_modes>1 
#==============================================================================
#   Example config.yaml
#==============================================================================
# ---
# name: set_0
# out: set_0.pdbqt
# in: conf_0.txt
# cpu: 4
# exhaustiveness: 12
# size:
# - 20
# - 20
# - 20
# rerun: 0
# overwrite_json: 1
# be_cutoff: -8.0
# dist_cutoff: 4.0
# scratch:  runs/ZINC_drugs_now_80
# in_json:  runs/ZINC_drugs_now_80/set_0.json
# out_json: runs/ZINC_drugs_now_80/set_0.json
# receptors:
# - /home/some/path/receptors/rec1.pdbqt
# - /home/some/path/receptors/rec2.pdbqt
# - /home/some/path/receptors/rec3.pdbqt
# - /home/some/path/receptors/rec4.pdbqt
# centers:
# - - -10.95
#   - 6.375
#   - 3.2
# - - 14.4
#   - -5.4
#   - 4.1
#
#==============================================================================
#   you can edit by hand! you can write it out with whatever scripting language you
#   are confortable with (ruby, python, perl all have yaml reading capabilities). In
#   perl, it's hashes and arrays.
#
#   name           => not used here, this used for pbs submission script.. 
#                     (-> pbs output for STDERR and STDOUT $name.o1790983 $name.e1790983)  
#   in             => name of the local config file that vina will load. 
#   out            => name of the local file that vina will dump configuration to. 
#                     This file is loaded and some info saved if be_cutoff satisfied.
#   rerun          => 0,1 will rerun if receptor and center keys already present (save if BE improves)
#   overwrite_json => 0,1 will die if trying to write to json file that exists
#   cpu            => integer, vina param number of cpus used in calculation 
#                     (converted to ppn for PBS submission)
#   exhaustiveness => integer, vina param how hard to search 
#   size           => how big of a box to use
#   be_cutoff      => number ("kcal/mol" vina), if BE < be_cutoff, store more info
#   dist_cutoff    => number (angstroms), if distance(atom_rec atom_lig) < dist_cutoff bin the 
#                     recepter residue, stored under NeighRes => {} 
#   scratch        => directory where the work is carried out 
#                     (may be in the same or different directory as json file)
#   in_json        => json containing ligands to run docking
#   out_json       => json file to write result from docking into info containing ligands
#   receptors      => array_ref [receptors to dock the ligands into at the centers]
#   centers        => array_ref of array_ref [ [x1,y1,z1], [x2,y2,z2] ]
#
#   This script will die unless the JSON file exists!
#
# OUTPUT:
#   STDERR: some carping if ligands are ignored or results overwritten
#   STDOUT: dumped yaml with info about best ligand, center, and receptor combo that was screened
#
#   if using pbs, this information will dump into the .e and .o files (.o accumulator script included)
#
#   JSON: This script stores a new json file with name based on the one provided. input.json -> input.json.new  
#   See merge_docks.pl for combining results.
#
#   See dock_these.pl for screening ligands/centers/receptors
#
use Modern::Perl;
use HackaMol;
use HackaMol::X::Vina;
use Math::Vector::Real;
use YAML::XS qw(Dump LoadFile);
use Time::HiRes qw(time);
use File::Slurp;
use JSON::XS;

my $tp1         = time;
my $tdock       = 0;
my $yaml_config = shift or die "pass yaml configuration file";
my $djob        = LoadFile($yaml_config);

#set some default configurations
#
$djob->{overwrite_json}=0           unless( exists( $djob->{overwrite_json}) );
$djob->{out_json}= $djob->{in_json} unless( exists( $djob->{out_json}) );
$djob->{rerun} = 0                  unless( exists( $djob->{rerun}) );

unless ($djob->{overwrite_json}){
  if( $djob->{out_json} eq $djob->{in_json}){
    die "trying to overwrite in_json, but overwrite_json set to 0";
  }
  if(-e $djob->{out_json}){
    die "trying to overwrite out_json, but overwrite_json set to 0";
  }
}

my $json_fn = $djob->{in_json};    #
die "$json_fn does not exist." unless ( -e $json_fn );
print STDERR "ligands are present, but will be ignored\n" if (exists($djob->{ligands}));

#set up vina docking object
my $vina = HackaMol::X::Vina->new(
    ligand         => 'tmp.pdbqt',
    receptor       => 'tmp.pdbqt',
    in_fn          => $djob->{in},
    out_fn         => $djob->{out},
    size           => V( @{$djob->{size}} ),
    cpu            => $djob->{cpu},
    exhaustiveness => $djob->{exhaustiveness},
    scratch        => $djob->{scratch},
    exe            => '~/bin/vina',
);


#load up the json file and ready the loop over existing ligands 
my $text = read_file( $json_fn, { binmode => ':raw' } );
my $json = new JSON::XS;
$json->incr_parse($text);

# hackamol instance for building and logging
my $hack = HackaMol->new( hush_read => 1, 
                             log_fn => Path::Tiny->tempfile(
                                         TEMPLATE => $djob->{name}."_XXXX",
                                         DIR      => $vina->scratch, 
                                         SUFFIX   => '.json',
                                         UNLINK   => 0,
                                        ),
                        );
my $ligand;
my $stor = {};

my $best = { BE => 0 };
my $fh = $hack->log_fn->openw_raw; 
my $doing_again = 0;
my $from_last_time = {};

while ( my $stor = $json->incr_parse ) {

    my ($lbase) = keys(%$stor);            # only 1 key
    my $ligand = $stor->{$lbase}{lpath};
    $vina->ligand($ligand);

    foreach my $receptor ( @{ $djob->{receptors} } ) {

        $vina->receptor($receptor);
        my $rec   = $hack->read_file_mol($receptor);
        my $rbase = $vina->receptor->basename('.pdbqt');

        $stor->{$lbase}{receptor}{$rbase}{rpath} = $vina->receptor->stringify;

        foreach my $center ( @{ $djob->{centers} } ) {

            $vina->center( V(@$center) );    # coercion needed here!!
            my $center_key = join( '_', @$center );
            if (
                exists( $stor->{$lbase}{receptor}{$rbase}{center}{$center_key} ) )
            {
                print STDERR "already docked $lbase $rbase $center_key\n";
                next if ($djob->{rerun} == 0);
                $doing_again = 1;
                print STDERR "overwriting $lbase $rbase $center_key if BE improves\n";
                $from_last_time = $stor->{$lbase}{receptor}{$rbase}{center}{$center_key};
            }

            my $t1  = time;
            my $mol = $vina->dock_mol(1);
            my $t2  = time;
            $tdock += $t2 - $t1;
            my $results =
              pack_up( $rec, $mol, $djob->{be_cutoff}, $djob->{dist_cutoff} );

            # replace new with old if new are not improvement
            if ($doing_again){
              if($results->{BE} > $from_last_time->{BE}){
                print STDERR "redocking did not improve\n";
                print STDERR "  last dock: ". $from_last_time->{BE}. " this dock: " . $results->{BE}."\n";
                $results = $from_last_time;
              }
              $doing_again = 0;
            }

            $stor->{$lbase}{receptor}{$rbase}{center}{$center_key} = $results;
            #best for ligand
            if ( $stor->{$lbase}{BEST}{BE} > $results->{BE} ) {
                $stor->{$lbase}{BEST}{BE}       = $results->{BE};
                $stor->{$lbase}{BEST}{COM}      = $results->{COM};
                $stor->{$lbase}{BEST}{receptor} = $rbase;
                $stor->{$lbase}{BEST}{center}   = $center_key;
            }
        }
    }
    #best overall
    if ( $best->{BE} > $stor->{$lbase}{BEST}{BE} ) {
         $best->{BE}       = $stor->{$lbase}{BEST}{BE};
         $best->{COM}      = $stor->{$lbase}{BEST}{COM};
         $best->{ligand}   = $lbase;
         $best->{TMass}    = $stor->{$lbase}{TMass};
         $best->{formula}  = $stor->{$lbase}{formula};
         $best->{receptor} = $stor->{$lbase}{BEST}{receptor};
         $best->{center}   = $stor->{$lbase}{BEST}{center};
    }
    print $fh encode_json $stor;
    #$hack->log_fn->append( encode_json $stor );
}

my $tp2 = time;
$best->{message}  = "results here are for the best observed overall for this set";
$best->{total_time} = $tp2 - $tp1;
$best->{dock_time}  = $tdock;

print Dump $best;

$hack->log_fn->move( $vina->scratch->child($djob->{out_json})->stringify );

sub pack_up {

    # save all that good stuff.  save more based on cutoff
    my ( $rec, $mlig, $be_cut, $dist_cut ) = ( shift, shift, shift, shift );
    $be_cut   = 0.0 unless defined($be_cut);
    $dist_cut = 4.0 unless defined($dist_cut);

    # results always returned
    my $results = {
        BE  => $mlig->get_score(0),
        COM => [ @{ $mlig->COM } ],
    };

    if ( $results->{BE} <= $be_cut ) {
        my @Zxyz = map { [ $_->Z, @{ $_->xyz } ] } $mlig->all_atoms;
        $results->{Zxyz} = \@Zxyz;
        my %Res;

        # accumulate binding site info
        foreach my $lat ( $mlig->all_atoms ) {

            my @pat = grep { $lat->distance($_) <= $dist_cut }
              grep { $_->Z != 1 } $rec->all_atoms;

            $Res{ $_->resname . "_" . $_->resid }++ foreach @pat;
        }
        $results->{NeighRes} = \%Res;

    }
    return $results;

}

