use Modern::Perl;
use PBS::Client;
use Path::Tiny;
use YAML::XS qw(DumpFile LoadFile Load Dump);


my $yaml = Load(<<'...');
---
cpu: 4 
exhaustiveness: 12
size:
- 20
- 20
- 20
overwrite_json: 1
rerun: 1
be_cutoff: -8.0
dist_cutoff: 4.0
scratch:  ligands/NCI_diversitySet2
receptors:
- receptors/some.pdbqt
centers:
- - -11.18
  - 0.06
  - -0.28
- - 9.85
  - 3.7
  - -5.55
- - 7.7
  - 3.2
  - 6.9
...

#adjust here if not all jsons need to be scanned
my @jsons = path($yaml->{scratch})->children(qr/\.json/);

foreach my $json (sort @jsons){
  $yaml->{name}     = $json->basename(qr/\.json/); 
  $yaml->{out}      = $yaml->{name} . ".pdbqt";
  $yaml->{in}       = $yaml->{name} . ".txt";
  $yaml->{in_json}  = $json->stringify;
  $yaml->{out_json} = $json->stringify;
  my $fyaml = path($yaml->{scratch})->child($yaml->{name} . ".yaml");
  DumpFile($fyaml,$yaml);

  my $client = PBS::Client->new();
  my $job    = PBS::Client::Job->new (
              queue => 'batch',
              name  => $yaml->{name},
              ppn   => $yaml->{cpu},
              nodes => 1,
              cput  => '24:00:00',
              cmd   => ["perl ligands_dock.pl $fyaml"],
  );
  say "perl ligands_dock.pl $fyaml";
  $client->qsub($job);
}

