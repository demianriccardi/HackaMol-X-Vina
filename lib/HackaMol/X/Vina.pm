package HackaMol::X::Vina;

#ABSTRACT: HackaMol extension for running Autodock Vina
use Moose;
use MooseX::StrictConstructor;
use Moose::Util::TypeConstraints;
use Math::Vector::Real;
use MooseX::Types::Path::Tiny qw(AbsPath) ;
use HackaMol; # for building molecules
use File::chdir;
use namespace::autoclean;
use Carp;

with qw(HackaMol::X::ExtensionRole);

has $_ => ( 
            is        => 'rw', 
            isa       => AbsPath, 
            predicate => "has_$_",
            required  => 1,
            coerce    => 1,
          ) foreach ( qw( receptor ligand ) );

has 'save_mol' => (
            is      => 'rw',
            isa     => 'Bool',
            default => 0,
);        
   

has $_ => (
    is        => 'rw',
    isa       => 'Num',
    predicate => "has_$_",
) foreach qw(center_x center_y center_z size_x size_y size_z);

has 'num_modes' => (
    is        => 'rw',
    isa       => 'Int',
    predicate => "has_num_modes",
    default   => 1,
    lazy      => 1,
);

has $_ => (
    is        => 'rw',
    isa       => 'Int',
    predicate => "has_$_",
) foreach qw(energy_range exhaustiveness seed cpu);


has 'center' => (
    is        => 'rw',
    isa       => 'Math::Vector::Real',
    predicate => "has_center",
    trigger   => \&_set_center,
);

has 'size' => (
    is        => 'rw',
    isa       => 'Math::Vector::Real',
    predicate => "has_size",
    trigger   => \&_set_size,
);

sub BUILD {
    my $self = shift;

    if ( $self->has_scratch ) {
        $self->scratch->mkpath unless ( $self->scratch->exists );
    }

    # build in some defaults
    $self->in_fn("conf.txt") unless ($self->has_in_fn);
    $self->exe("~/bin/vina") unless $self->has_exe;

    unless ( $self->has_out_fn ) {
      my $outlig = $self->ligand->basename;
      $outlig =~ s/\.pdbqt/\_out\.pdbqt/;
      $self->out_fn($outlig); 
    }

    unless ( $self->has_command ) {
        my $cmd = $self->build_command;
        $self->command($cmd);
    }

    return;
}

sub _set_center {
    my ( $self, $center, $old_center ) = @_;
    $self->center_x( $center->[0] );
    $self->center_y( $center->[1] );
    $self->center_z( $center->[2] );
}

sub _set_size {
    my ( $self, $size, $old_size ) = @_;
    $self->size_x( $size->[0] );
    $self->size_y( $size->[1] );
    $self->size_z( $size->[2] );
}

#required methods
sub build_command {
    my $self = shift;
    my $cmd;
    $cmd = $self->exe;
    $cmd .= " --config " . $self->in_fn->stringify;

    # we always capture output
    return $cmd;
}

sub _build_map_in {
    # this builds the default behavior, can be set anew via new
    return sub { return ( shift->write_input ) };
}

sub _build_map_out {
    # this builds the default behavior, can be set anew via new
    my $sub_cr = sub {
        my $self = shift;
        my $qr   = qr/^\s+\d+\s+(-*\d+\.\d)/;
        my ( $stdout, $sterr ) = $self->capture_sys_command;
        my @be = map { m/$qr/; $1 }
          grep { m/$qr/ }
          split( "\n", $stdout );
        return (@be);
    };
    return $sub_cr;
}

sub dock {
  my $self      = shift;
  my $num_modes = shift;
  $self->num_modes($num_modes) if defined($num_modes);
  $self->map_input;
  return $self->map_output;
}

sub dock_mol {
  # want this to return configurations of the molecule
  my $self      = shift;
  my $num_modes = shift;
  $self->num_modes($num_modes) if defined($num_modes);
  $self->map_input; 
  local $CWD = $self->scratch if ( $self->has_scratch );
  my @bes = $self->map_output; # this is fragile... broken if map_out changed...
  my $mol = HackaMol -> new(hush_read => 1)
                     -> read_file_mol($self->out_fn->stringify);
  $mol->push_score(@bes);
  return ($mol);
}

sub write_input {
    my $self = shift;
    my $input;
    $input .= sprintf( "%-15s = %-55s\n", 'out', $self->out_fn->stringify );
    $input .= sprintf( "%-15s = %-55s\n", 'log', $self->log_fn->stringify )
      if $self->has_log_fn;
    foreach my $cond (
        qw(receptor ligand cpu num_modes energy_range exhaustiveness seed))
    {
        my $condition = "has_$cond";
        $input .= sprintf( "%-15s = %-55s\n", $cond, $self->$cond )
          if $self->$condition;
    }
    foreach my $metric (qw(center_x center_y center_z size_x size_y size_z)) {
        $input .= sprintf( "%-15s = %-55s\n", $metric, $self->$metric );
    }
    $self->in_fn->spew($input);
    return ($input);
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 SYNOPSIS

  use Modern::Perl;
  use HackaMol;
  use HackaMol::X::Vina;
  use Math::Vector::Real;

  my $receptor = "receptor.pdbqt";
  my $ligand   = "lig.pdbqt",
  my $rmol     = HackaMol -> new( hush_read=>1 ) -> read_file_mol( $receptor );
  my $lmol     = HackaMol -> new( hush_read=>1 ) -> read_file_mol( $ligand );
  my $fh = $lmol->print_pdb("lig_out.pdb");

  my @centers = map  {$_ -> xyz}
                grep {$_ -> name    eq "OH" }
                grep {$_ -> resname eq "TYR"} $rmol -> all_atoms;

  foreach my $center ( @centers ){

      my $vina = HackaMol::X::Vina -> new(
          receptor       => $receptor,
          ligand         => $ligand,
          center         => $center,
          size           => V( 20, 20, 20 ),
          cpu            => 4,
          exhaustiveness => 12,
          exe            => '~/bin/vina',
          scratch        => 'tmp',
      );

      my $mol = $vina->dock_mol(3); # fill mol with 3 binding configurations

      printf ("Score: %6.1f\n", $mol->get_score($_) ) foreach (0 .. $mol->tmax); 

      $mol->print_pdb_ts([0 .. $mol->tmax], $fh); 

    }

    $_->segid("hgca") foreach $rmol->all_atoms; #for vmd rendering cartoons.. etc
    $rmol->print_pdb("receptor.pdb");


=head1 DESCRIPTION

HackaMol::X::Vina provides an interface to AutoDock Vina. This class does not include the AutoDock Vina program, which is 
<released under a very permissive Apache license|http://vina.scripps.edu/manual.html#license>, with few restrictions on 
commercial or non-commercial use, or on the derivative works, such is this. Follow these 
<instructions | http://vina.scripps.edu/manual.html#installation> to acquire the program. Most importantly, if you use this 
interface effectively, please be sure to cite AutoDock Vina in your work:

O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461 

Since HackaMol has no pdbqt writing capabilities (yet, HackaMol can read pdbqt files), the user is required to provide
those  files. This is still a work in progress and the API may still change. Documentation will improve as API
gets more stable... comments welcome!

=method write_input

This method takes no arguments; it returns, as a scalar, the input constructed from attributes.  This method is called by map_input method via the map_in attribute. 

=method map_input

provided by HackaMol::X::ExtensionRole. Writes the configuration file for Vina. See dock and dock_mol methods.

=method map_output

provided by HackaMol::X::ExtensionRole. By default, this method returns the docking scores as an array.

=method dock_mol

this method takes the number of binding modes (Integer) as an argument (Int). The argument is optional, and the num_modes attribute is rewritten if passed. This method calls the map_input and map_output methods for preparing and running Vina. It loads the resulting pdbqt and scores into a HackaMol::Molecule object.  The scores are stored into the score attribute provided by the HackaMol::QmMolRole. See the synopsis for an example.

=method dock

this method is similar to dock_mol, but returns only the scores.


