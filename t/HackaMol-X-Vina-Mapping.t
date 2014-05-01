#!/usr/bin/env perl

use strict;
use warnings;
use Test::Moose;
use Test::More;
use Test::Fatal qw(lives_ok dies_ok);
use Test::Dir;
use Test::Warn;
use HackaMol::X::Vina;
use HackaMol;
use Math::Vector::Real;
use File::chdir;
use Path::Tiny;
use Cwd;

my $cwd      = getcwd;
my $receptor = path('t/lib/receptor.pdbqt'); 
my $ligand   = path('t/lib/ligand.pdbqt' ); 
my $rmol     = HackaMol->new(hush_read=>1)->read_file_mol($receptor);    
my $lmol     = HackaMol->new(hush_read=>1)->read_file_mol($ligand);    
my $ligout   = $ligand->basename;
$ligout      =~ s/\.pdbqt/_new\.pdbqt/;       

my $obj = HackaMol::X::Vina->new(
        receptor       => $receptor->absolute->stringify,  # needs coercion!
        ligand         => $ligand->absolute->stringify,
        in_fn          => "conf.txt",
        out_fn         => $ligout,
        center         => V(6.865, 3.449, 85.230),
        size           => V(10,10,10),
        cpu            => 1,
        num_modes      => 5,
        exhaustiveness => 1,
        exe            => '~/bin/vina',
        scratch        => 't/tmp',
        seed           => 314159,
    );

my $input = $obj->map_input;
print $input;

my @bes   = $obj->map_output;
my @Be_expected = qw(
  -6.6
  -5.9
  -5.8
  -5.7
  -5.1
);

is_deeply(\@bes, \@Be_expected, 'binding energies computed with vina');

$obj->center( V( 18.073, -2.360, 90.288 ) );
$input = $obj->map_input;
print $input;
@bes   = $obj->map_output;
print $_ . "\n" foreach @bes;

done_testing();

