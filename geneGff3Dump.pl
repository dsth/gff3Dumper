#!/usr/local/bin/perl5.8.0

# Daniel S. T. Hughes
# dsth@ebi.ac.uk
# dsth@cpan.org
# dsth@cantab.net

# for local use comment out
#use lib '/home/dsth/ensembl_source/ensembl_56/modules/';

# perl geneGff3Dump.pl -dbuser ensro -dbhost mysql-eg-live-1.ebi.ac.uk -dbport 4159 -dbname anopheles_gambiae_core_5_58_3l

use strict;
use warnings;
no warnings 'uninitialized';
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::MetaContainer;
use Getopt::Long;
use Bio::SeqIO;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

# Variaveis ={{{1
my $dbhost = 'mysql-eg-live-1.ebi.ac.uk';
my $noutr = 0;
my $dbport = 4159;
my $dbuser = 'ensro';
my $dbname = 'anopheles_gambiae_core_4_58_3l';
my $metazoa = 0;
#y don't like ns collisions. file scoped so diff name and re-declare in lower scoping
my $source_db = q{};

my $molecule_type     = "dsDNA";
my $dbxref_name;

my $verbose;
my $verbose_xref;
my $test;
my $chrom_file;
my $data_dir = '.';
my $dbxref_text;
my $gff_header;
my $output_file; # = 'gene.gff3';
#=}}}

my %SO = ( 
# SO identifiers ={{{1
    'gene'       => 'SO:0000704',
    'exon'       => 'SO:0000147',     
    'CDS'        => 'SO:0000316',              
    '5UTR'       => 'SO:0000204',
    '3UTR'       => 'SO:0000205',
    'transcript' => 'SO:0000234'     # SO code for mRNA
#=}}}
);

#y DB synonyms - i.e. don't prinki full long EnsEMBL names for them
my %xref_database_aliases = (
#={{{1 xref aliases
    'WikiGene'                  => 'WikiGene',    'Uniprot/SPTREMBL_predicted'  => 'UniProtKB',
    'Uniprot/SPTREMBL'          => 'UniProtKB',   'Uniprot/SWISSPROT'           => 'UniProtKB',
    'protein_id_predicted'      => 'protein_id',  'protein_id'                  => 'protein_id',
    'EMBL_predicted'            => 'GenBank',     'EMBL'                        => 'GenBank',
    'UniGene'                   => 'UniGene',     'RefSeq_dna_predicted'        => 'RefSeq_NA',
    'RefSeq_peptide_predicted'  => 'RefSeq_Prot',
    # Aedes
    'AedesGenBank'              => 'protein_id',
#=}}}
);

# A lookup of an element in a hash is very fast, while grep and the likes have to go through the whole array, which is slower.
#y is hash table look up any faster than grep
my %xref_databases_ignored = (
# xrefs to ignore ={{{2
    'AFFY_Plasmodium_Anopheles'     => 1,  'ARRAY_EMBL_MMC2_12k_v1'        => 1,
    'ARRAY_EMBL_MMC1_20k_v1'        => 1,  'ARRAY_LIV_GAMDETOX_0_25k_v1'   => 1,
    'ARRAY_LIV_GAMDETOX_0_25k_v2'   => 1,  'ARRAY_JHSPH_AGGAMBER_15k_v1'   => 1,
    'Celera_Gene'                   => 1,  'Celera_Trans'                  => 1,
    'Celera_Pep'                    => 1,  'GO'                            => 1,
    'Uniprot/Varsplic'              => 1,  'PDB'                           => 1,
    'EntrezGene'                    => 1,
    # Anopheles
    'agam-detox-1J20'               => 1,  'agam-detox-1C18'               => 1,
    'agam-detox-1C22'               => 1,  'agam-detox-1C6'                => 1,
    'agam-detox-1C2'                => 1,  'agam-detox-1C3'                => 1,
    'agam-detox-1E1'                => 1,  'agam-detox-1M6'                => 1,
    'agam-detox-1B19'               => 1,  
    # Aedes
    'ARRAY_UCR_GillMgMT_0_2K_v2'    => 1,  'ARRAY_ND_TIGRTC_9_6K_v1'       => 1,
    'ARRAY_LIV_AEGDETOX_0_25k_v1'   => 1,  
#=}}}
);

&GetOptions (
# Options ={{{1
    #y option
    '-metazoa'       => \$metazoa,
    '-noutr'         => \$noutr,
    #y args
    #b duh forgot to make it arg string and not option with :s
    '-source:s'        => \$source_db,
    '-dbhost=s'      => \$dbhost,
    '-dbport=i'      => \$dbport,
    '-dbname=s'      => \$dbname,
    '-dbuser=s'      => \$dbuser,
    '-data_dir=s'    => \$data_dir,
    '-contig_file=s' => \$chrom_file,
    '-output=s'      => \$output_file,
    '-test'          => \$test,
    '-verbose'       => \$verbose,
    '-verbose_xref'  => \$verbose_xref,
#=}}}
);

Bio::EnsEMBL::Registry->load_registry_from_db(
-host   => $dbhost,
-port   => $dbport,
-user   => $dbuser,
-dbname => $dbname,
);

if (!$metazoa) { 
    $output_file = $dbname.'.gff3';
    &_single_db($dbhost, $dbport, $dbuser, $dbname, $output_file, $source_db);
}
else {

    print qq{\nChecking for Metazoan databases.\n\n};
    my $metazoans =  &_metazoa_list($dbhost, $dbport, $dbuser, $dbname);

    die qq{\nCould not find any metazoa databases.} if (!@{$metazoans});
    #die qq{\nCould not find any metazoa databases} if (scalar @{$metazoans} == 0);

    print 'Found ', scalar @{$metazoans}, qq{ Metazoan databases:\n};

    print qq{\n* Species: }, $_->{alias}, ' (', $_->{common_name}, 
      '), assembly: ', $_->{assembly}, ', source: ', $_->{source} for (@{$metazoans});

    print qq{\n};

    for my $species (@{$metazoans}) {

        my $alias = $species->{alias};

        print qq{\nChecking adaptor for \'$alias\' database\n};

        my $ad = $species->{adaptor};

        $ad->dbc->connect or die qq{\nAdaptor DB connection gone} if (!$ad->dbc->connected);

        #my $source = $slice->{source} ? $...
        my $source = $species->{source};
        $source ||= 'Ensembl Genomes';

        #$ad = Bio::EnsEMBL::Registry->get_DBAdaptor('Culex quinquefasciatus', 'Core') 
        #my $ad = Bio::EnsEMBL::Registry->get_adaptor($alias, 'Core', 'Slice') 
         # or die qq{\ncannot make adpator};

        my $file = $alias;

        $file =~ s/[\ ]+/_/g;

        $file .= '.gff3';

        print qq{\nOpening file \'$data_dir/$file\' for writting\n};
        open my $met_hndl, '>', qq{$data_dir/$file} or die qq{\nCannot open file $file for writting};
        
        print STDOUT "Writting to output file $data_dir/$file\n\n";

        my @slices = @{$ad->get_SliceAdaptor->fetch_all('toplevel')};
        #my @slices = @{$ad->fetch_all('toplevel')};

        die qq{\nCould not fetch slices} if (scalar @slices == 0);

        &_process_slice(\@slices, $met_hndl, $source);

        close $met_hndl;
    }
}




sub _single_db {
    
    my ($dbhost, $dbport, $dbuser, $dbname, $output_file, $source) = @_;

    # open database connection {{{1
    my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( 
        -host   => $dbhost,
        -port   => $dbport,
        -user   => $dbuser,
        -dbname => $dbname,
    );

    #b now that the two methods are converging this is pretty ineffcient

    $source ||= $db->get_MetaContainer->list_value_by_key('species.division')->[0];

    open my $OUT_hndl, '>', qq{$data_dir/$output_file};
    print STDOUT "// Output file $data_dir/$output_file\n\n";

    #y don't need SliceAd persistence... 
    my $slice_ad = $db->get_SliceAdaptor() or die qq{\nCould not create slice apadtor};

    my @slices = @{$slice_ad->fetch_all('toplevel')};
    die qq{\nCould not fetch slices} if (scalar @slices == 0);

    &_process_slice(\@slices, $OUT_hndl, $source);

    close $OUT_hndl;

    return;
}

sub _metazoa_list {

    my ($dbhost, $dbport, $dbuser, $dbname) = @_;

    #y/ must redeclare any file-scoped var you don't want invdaing - or use diff name

    my @metazoa;

    for my $ad (@{Bio::EnsEMBL::Registry->get_all_DBAdaptors}) {

        my $meta = $ad->get_MetaContainer;

        if ($meta->list_value_by_key('species.division')->[0] eq 'EnsemblMetazoa') {

            push @metazoa, +{
                #y just take the first alias and use it form DB construction
                alias => $meta->list_value_by_key('species.alias')->[0],
                schema => $meta->list_value_by_key('schema_version')->[0],
                assembly => $meta->list_value_by_key('assembly.default')->[0],
                common_name => $meta->list_value_by_key('species.common_name')->[0],
                source => $meta->list_value_by_key('provider.name')->[0],
                adaptor => $ad,
                #$meta->list_value_by_key('species.ensembl_alias_name')->[0],
                #print qq{\nstable id prefix: }, $meta->list_value_by_key('species.stable_id_prefix')->[0];
                #print qq{\nclassification - multivalue: }, @{$meta->list_value_by_key('species.classification')};
            #print qq{\nalias - multivalue: }, @{$meta->list_value_by_key('species.alias')};
            };
    
        }
        else {
            #/ twat this was killing the program?!?
            #y we can retive the DBI object handle dicrecly of through connection
            $ad->db_handle->disconnect();
            #$ad->dbc->db_handle->disconnect();
        }
    }

    return \@metazoa;
}

sub _process_slice {

    my ($slices, $OUT_hndl, $source) = @_;

    SEQIDLOOP:
    for my $i (0..$#{$slices}) { 

        my $slice = $slices->[$i];

        my $name = $slice->name;
        my $seqid = $slice->seq_region_name;
        my $acc = $slice->accession_number;
        my $type = $slice->coord_system_name;
        my $start = $slice->start();
        my $end   = $slice->end();

        #y print header for file
        print $OUT_hndl $i == 0 
          ? &_theredoc(\%SO, $dbname, $name, $source) 
          : "##sequence-region $seqid\n";

        print STDOUT qq{Slice: $seqid\n};

        #y feature entry for seqid/contig - whatever toplevel is in this case
        print $OUT_hndl qq{$name\t$source\tcontig\t$start\t$end\t.\t.\t.\t}
          #b attributes column
          . qq{ID=$seqid;}
          . q{molecule_type=dsDNA;}  
          . qq{GenBank:$acc;}
          . q{translation_table=1;} 
          . q{topology=linear;}
          . qq{localization=chromosomal;\n}; 
        
        #g call gene processing sub
        &_process_gene($OUT_hndl, $slice, $seqid, $source);

    } 
   return;
}

sub _process_gene {
    
    my ($OUT_hndl, $slice, $seqid, $source) = @_;

    for my $gene (@{$slice->get_all_Genes}) { 

#y/ temp
#next if ($gene->stable_id ne 'AGAP000201');

        print STDOUT "\n// Gene ",$gene->stable_id,"\t" if $verbose;

        #y set $strand for whole gene - its in ensembl it should always agree...
        my $strand = $gene->strand == 1 ? '+' : '-'; 

        print STDOUT ' Gene ', $gene->stable_id, qq{\n};

        #y print the initial feature for the gene
        print $OUT_hndl qq{$seqid\t$source\tgene\t}, $gene->start,qq{\t}, $gene->end, 
          #qq{\t.\t$strand\t.\tID=}, $gene->stable_id, qq{;web_id=}, $gene->stable_id, qq{;\n};
          qq{\t.\t$strand\t.\tID=}, $gene->stable_id, qq{;\n};

        #y process gene's transcripts
        TRANSLOOP: # {{{1
        for my $trans (@{$gene->get_all_Transcripts}) {
            
            #b bin contaminants
            next TRANSLOOP if ($trans->biotype() eq "bacterial_contaminant");

            print STDOUT qq{  Transcript }, $trans->stable_id, qq{\n} if $verbose;

            #y we are using protein_coding and NOT translatable as condition to avoid lossing crappy things but print message
            my $stuff = $trans->biotype eq "protein_coding" 
                ? ['mRNA','stable_id='.$trans->stable_id.';']
                : [$trans->biotype, q{}];

            #y print transcript feature cols 1-8
            print $OUT_hndl qq{$seqid\t$source\t}, $stuff->[0], qq{\t}, $trans->start, qq{\t}, $trans->end, qq{\t.\t$strand\t.\t};

            #y print transcript attributes col (minus xrefs and descriptions)
            print $OUT_hndl qq{ID=}, $trans->stable_id, qq{;Parent=}, $gene->stable_id, qq{;}; #, $stuff->[1];
                #qq{;web_id=}, $trans->stable_id, 

            #r for () =/= for (0)

            #g get and print xrefs using _process_xrefs - if none, clearly prints empty string
            print $OUT_hndl &_process_xref($OUT_hndl, $trans);
            #print $OUT_hndl $xref_string.';' if ($xref_count != 0);

            #y print description
            print $OUT_hndl $trans->description ? 'description='.$trans->description.qq{\n} : qq{description=hypothetical protein;\n};

            #g call process exons to print exons and when protein_coding CDS too        
            &_process_exons($OUT_hndl, $seqid, $gene, $trans, $strand, $source);
            
        } 
    }
   return;
}

sub _process_exons {

    my ($OUT_hndl, $seqid, $gene, $trans, $strand, $source) = @_;

    my $cds_count = 0;
    my $phase_minusOne;   
    my $length_minusOne;

    my $exon_string = q{};
    my $cds_string = q{};
    my $fiveprime = qq{};
    my $threeprime = qq{};

    #y set translation flag - i.e. it has CDS processing and get translation object
    my $transl = $trans->translation ? $trans->translation : 0; # 
    my $ending = 0;
    #my $transl = $trans->biotype eq 'protein_coding' ? $trans->translation : 0; # 

    #y define scalars to store important anchors
    my $first_transl_exon;
    my $last_transl_exon; 
    my $transl_start; 
    my $transl_end;
    my $physical_start;
    my $physical_end;
    my $utr_5_exists = 1;
    my $utr_3_exists = 1;
            

    #y if its protein_coding set these vars for later reference
    if ($transl) {

        #$first_exon  = $trans->start_Exon;
        #$last_exon    = $trans->end_Exon;
        $physical_start = $strand eq '+' ? $trans->start : $trans->end;
        $physical_end = $strand eq '+' ? $trans->end : $trans->start;

        $first_transl_exon  = $transl->start_Exon;
        $last_transl_exon    = $transl->end_Exon;
        #r these are relative to the exon start/end
        $transl_start     = $transl->start;
        $transl_end       = $transl->end;

    }

    #y iterate through the exons
    for my $exon (@{$trans->get_all_Exons})  {

        #r no need for ordering etc., as get_all_Exons returns all exons in the order of their translation - i.e. same irrespective of strand
        # for my $exon (@ordered_exons) { 

        #r in list matching returns ALL the match vars $1, $2... this can be done by assignment...
        #my $isoform = $trans->stable_id  =~ (/\S+\-R(\S+)/);
        my $isoform = q{};
        if ($trans->stable_id  =~ (/\S+\-R(\S+)/)) { $isoform = $1 }

        #y get the basic info
        my @exon_position = ($exon->start, $exon->end);

        $exon_string .= qq{$seqid\t$source\texon\t$exon_position[0]\t$exon_position[1]\t.\t$strand\t.\tID=}.
        $exon->stable_id.qq{$isoform;Parent=}.$trans->stable_id.qq{;\n};
        #$exon_string .= qq{$seqid\tVectorBase\texon\t$exon_start\t$exon_end\t.\t$strand\t.\tID=}.

        #y CDS handling
        if ($transl) {

            my  $alt_strand = $gene->strand == 1 ? '-' : '+';
            my  $switch = $gene->strand == 1 ? 0 : 1 ;

            #r these are objects - i.e. in condition will be evaulated by memory location hence use string
            #r could just use phase == -1 etc to distinguish...

            #b silly utr crap so need to print five prime...
            #f ordered list so if its not the first and there's nothing in the string we're still 5' to start site
            #next if ($exon ne $first_transl_exon && !$cds_string);
            if ($exon ne $first_transl_exon && !$cds_string) {

                #y literally just concat the exon info again?!? but this time call it utr
                $fiveprime .= qq{$seqid\t$source\tfive_prime_utr\t$exon_position[0]\t$exon_position[1]\t.\t$strand\t.\tParent=}
                  . $trans->stable_id.qq{;\n};
                  #y/ done for this exon
                next;
            }


            #y override start/end for first/last exons
            #r if -1 adjust end of first and start of last
            #r if +1 adjust start of first and end of last

            my $l = $exon->length;

            #/ need to have a check that the original value doesn't equal old in which case there's no utr and switch flag

            #y now doing utrs as well
            #$exon_position[$switch] = eval(qq{$exon_position[$switch]$strand$transl_start${alt_strand}1}) if ($exon eq $first_transl_exon);
            if ($exon eq $first_transl_exon) {

                my $hack = $exon_position[$switch];

                my @utr_hack;

                #b/ dumb ass utr is literally just the piece between original location and new
                $utr_hack[$switch] = $exon_position[$switch];

                #y/ the actual code
                #y modify start and end positions for first and last coding exons for the CDS
                #b this interpolates strand +/- into string and uses eval to force its execution as code (force string as operator).
                $exon_position[$switch] = eval(qq{$exon_position[$switch]$strand$transl_start${alt_strand}1});

                # this should be if else - i.e. if there's no utr no need for calc
                #b annoying me so just use whether the value changes or not
                #$utr_5_exists = 0 if ($exon_position[1-$switch] == $physical_start); #it should equal start w/o utr?!?
                $utr_5_exists = 0 if ($hack == $exon_position[$switch]);

                #b hack for five prime utr start - eew.
                #print qq{\ntest: }, $trans->start;
                #$fiveprime .= $strand eq q{+} ? $trans->start.qq{\t} : $trans->end.qq{\t};
                #$fiveprime .= $strand eq q{+} ? $trans->start.qq{\t} : $trans->end.qq{\t};
                #$utr_hack[$switch] = $physical_start;
                #duh: $utr_hack[1-$switch] = $exon_position[1-$switch]; # i.e. one stays the same

                #r using +/- string as operator
                $utr_hack[1-$switch] = eval(qq{$exon_position[$switch]${alt_strand}1});

                $fiveprime .= qq{$seqid\t$source\tfive_prime_utr\t$utr_hack[0]\t}
                  .qq{$utr_hack[1]\t.\t$strand\t.\tParent=}.$trans->stable_id.qq{\n};

            }

            if ($exon eq $last_transl_exon) {

                my $hack = $exon_position[1-$switch];

                my @utr_hack;

                #b/ dumb ass utr is literally just the piece between original location and new
                $utr_hack[1-$switch] = $exon_position[1-$switch];

                #y/ actual code
                $exon_position[1-$switch] = eval(qq{$exon_position[1-$switch]$strand($transl_end-$l)});

                $utr_3_exists = 0 if ($hack == $exon_position[1-$switch]);
                #$utr_3_exists = 0 if ($exon_position[$switch] == $physical_end);

                $utr_hack[$switch] = eval(qq{$exon_position[1-$switch]${strand}1});

                $threeprime .= qq{$seqid\t$source\tthree_prime_utr\t$utr_hack[0]\t}
                  .qq{$utr_hack[1]\t.\t$strand\t.\tParent=}.$trans->stable_id.qq{\n};

            }

            #$exon_position[1-$switch] = eval(qq{$exon_position[1-$switch]$strand($transl_end-$l)}) if ($exon eq $last_transl_exon);

            #y adjust length for first and last CDS
            $l = $exon_position[1] - $exon_position[0] + 1 if ($exon eq $first_transl_exon || $exon eq $last_transl_exon);

            #f these are ORDERED EXONS!?! so its always the same irrespective of strand
            #my $p = $exon->phase == -1 ? 0 : (3 - $exon->phase) % 3;
            my $p = $cds_count == 0 ? 0 : ((3 -(($length_minusOne%3)-$phase_minusOne) ) % 3);

            $cds_string .= qq{$seqid\t$source\tCDS\t$exon_position[0]\t$exon_position[1]\t.\t$strand\t$p\t}
                #.qq{Parent=}.$trans->stable_id.qq{;\n}; 
                .qq{Parent=}.$trans->stable_id.qq{;\n} if (!$ending); 

            #y stash values for phase calculations for next iteration
            $length_minusOne = $l;
            $phase_minusOne = $p;
            $cds_count++;

            #b/ over-ride this for pointless utr handling - instead flick ending switch
            # stop entering CDS handling after the last coding exon
            #$transl = 0 if ($exon eq $last_transl_exon);
            #y instead print if last exon already past using below switch
            $threeprime .= qq{$seqid\t$source\tthree_prime_utr\t$exon_position[0]\t$exon_position[1]\t.\t$strand\t.\tParent=}
              . $trans->stable_id.qq{;\n} if ($ending);

#y/ stop printing CDS now
            $ending = 1 if ($exon eq $last_transl_exon);

        } #else {} # phase is -1 } but CDS have no phase

    }

    print $OUT_hndl $exon_string;

    #b now compound condition as no point in printing if doesn't exist
    print $OUT_hndl $fiveprime if (!$noutr && $utr_5_exists);

    print $OUT_hndl $cds_string;

    #y being naughty-ish and using the file scoped lexical directly
    #print $OUT_hndl $threeprime if !$noutr;
    #y start using short-circuits more often
    #$noutr or print $OUT_hndl $fiveprime.qq{\t.\t$strand\t.\t\n};
    #$noutr or print $OUT_hndl $threeprime.qq{\t.\t$strand\t.\t\n};
    
    print $OUT_hndl $threeprime if (!$noutr && $utr_3_exists);

    return
}

sub _process_xref {
    my ($OUT_hndl,$trans) = @_;

    my $xref_string = q{Dbxref=};
    my $xref_count = 0;
    my $gene_symbol;

    for my $xref (@{$trans->get_all_DBLinks}) {

        $dbxref_name = $xref->dbname;
        #y ignore blacklisted DBs
        next if ($xref_databases_ignored{$dbxref_name}); # ignore some DBs

        #y no longer skip unknown dbs whitelist is now just alias list
        # skipping unknown DBs
        # next if (!$xref_database_aliases{$dbxref_name});
        #y instead if the alias exists that's its name else its the original name
        $dbxref_name = exists $xref_database_aliases{$dbxref_name}
          ? $xref_database_aliases{$dbxref_name}
          : $dbxref_name;

        #y anopheles has gene symbol stored as xref
        # we haven't actually printed the xrefs yet so can do this in middle
        if ($dbxref_name eq 'Anopheles_symbol' ) {
            print $OUT_hndl 'gene_symbol=', $xref->primary_id, ';';
            next;
        }

        $xref_string .= ',' if ($xref_count != 0); # comma before next pair
        $xref_string .= $dbxref_name . ':' . $xref->primary_id;
        $xref_count++;
    }
            
    return $xref_string.';' if ($xref_count != 0);
    return;
}

sub _theredoc {
my ($SO, $dbname, $seqid, $start, $end, $source) = @_;
my %SO = %{$SO};
return <<"HEADER";
##gff-version 3
##feature-ontology so.obo
##attribute-ontology gff3_attributes.obo
#
# Dumped from $source database.
# 
# using $SO{gene} for VectorBase Gene
# using $SO{transcript} for VectorBase Transcript
# using $SO{exon} for VectorBase Exon
# using $SO{CDS} for VectorBase CDS
# using $SO{'5UTR'} for VectorBase 5'UTR 
# using $SO{'3UTR'} for VectorBase 3'UTR
# 
##sequence-region $seqid
HEADER
}

