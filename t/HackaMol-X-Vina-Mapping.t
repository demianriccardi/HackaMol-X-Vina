#!/usr/bin/env perl

use strict;
use warnings;
use Test::More;
use Test::Dir;
use HackaMol::X::Vina;
use HackaMol;
use Math::Vector::Real;
use Path::Tiny;

my $receptor = path('t/lib/receptor.pdbqt'); 
my $ligand   = path('t/lib/ligand.pdbqt' ); 
my $rmol     = HackaMol->new(hush_read=>1)->read_file_mol($receptor);    
my $lmol     = HackaMol->new(hush_read=>1)->read_file_mol($ligand);    
my $ligout   = $ligand->basename;
$ligout      =~ s/\.pdbqt/_out\.pdbqt/;       

my $obj = HackaMol::X::Vina->new(
        receptor       => $receptor->absolute->stringify, 
        ligand         => $ligand->absolute->stringify,
        in_fn          => "conf.txt",
        out_fn         => $ligout,
        center         => V(  6.865, 3.449, 85.230 ),
        size           => V( 10, 10, 10 ),
        cpu            => 1,
        num_modes      => 2,
        exhaustiveness => 1,
        exe            => '~/bin/vina',
        scratch        => 't/tmp',
        seed           => 314159,
);

my $input = $obj->map_input;
my @bes   = $obj->map_output;
my @Be_expected = qw(
  -6.6
  -5.9
);

is_deeply(\@bes, \@Be_expected, 'binding energies computed with vina');

$obj->center( V( 18.073, -2.360, 90.288 ) );

$input = $obj->map_input;
@bes   = $obj->map_output;
@Be_expected = qw(
-4.2
-3.9
);

is_deeply(\@bes, \@Be_expected, 'binding energies computed with vina, new center');

$obj->scratch->remove_tree;
dir_not_exists_ok( "t/tmp", 'scratch directory deleted' );

done_testing();

