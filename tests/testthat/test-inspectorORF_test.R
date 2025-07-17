test_that("importing data",
{
  # ms_data_file <- "../../../Annotation/control_MS_evidence.bed"
  #
  # bed_test_2 <- "../../../Annotation/test2.bed"
  # bed_test <- "../../../Annotation/test.bed"
  #
  # gtf_test <- "../../../Annotation/gencode.v30.annotation.gtf"
  # gtf_test_2 <- "../../../Annotation/gencode.v44.annotation.gtf"
  #
  # genome_test <- "../../../Annotation/GRCh38.p12.primary_assembly.genome.2bit"
  # genome_test_2 <- "../../../Annotation/GRCh38.p14.primary_assembly.genome.2bit"
  #
  # test_tracks <- inspectorORF::create_tracks(bed_test_2, gtf_test_2, genome_test_2, ms_data_file)
})

test_that("finding ORFs",
{
  ms_data_file <- "../../../Annotation/control_MS_evidence.bed"

  bed_test_2 <- "../../../Annotation/test2.bed"
  bed_test <- "../../../Annotation/test.bed"

  gtf_test <- "../../../Annotation/gencode.v30.annotation.gtf"
  gtf_test_2 <- "../../../Annotation/gencode.v44.annotation.gtf"

  genome_test <- "../../../Annotation/GRCh38.p12.primary_assembly.genome.2bit"
  genome_test_2 <- "../../../Annotation/GRCh38.p14.primary_assembly.genome.2bit"

  test_tracks <- inspectorORF::create_tracks(bed_test_2, gtf_test_2, genome_test_2, ms_data_file)

  test_orfs <- inspectorORF::find_orfs(test_tracks, "ENST00000643391.1")
  test_orfs
})

test_that("plotting transcript",
{
  # ms_data_file <- "../../../Annotation/control_MS_evidence.bed"
  #
  # bed_test_2 <- "../../../Annotation/test2.bed"
  # bed_test <- "../../../Annotation/test.bed"
  #
  # gtf_test <- "../../../Annotation/gencode.v30.annotation.gtf"
  # gtf_test_2 <- "../../../Annotation/gencode.v44.annotation.gtf"
  #
  # genome_test <- "../../../Annotation/GRCh38.p12.primary_assembly.genome.2bit"
  # genome_test_2 <- "../../../Annotation/GRCh38.p14.primary_assembly.genome.2bit"
  #
  # test_tracks <- inspectorORF::create_tracks(bed_test_2, gtf_test_2, genome_test_2, ms_data_file)
  #
  # test_transcript_plot <- inspectorORF::transcript_plot(test_tracks, "ENST00000486256.5", with_ORF_start_position = 1702)
  # test_transcript_plot
})

test_that("plotting ORF",
{
  # ms_data_file <- "../../../Annotation/control_MS_evidence.bed"
  #
  # bed_test_2 <- "../../../Annotation/test2.bed"
  # bed_test <- "../../../Annotation/test.bed"
  #
  # gtf_test <- "../../../Annotation/gencode.v30.annotation.gtf"
  # gtf_test_2 <- "../../../Annotation/gencode.v44.annotation.gtf"
  #
  # genome_test <- "../../../Annotation/GRCh38.p12.primary_assembly.genome.2bit"
  # genome_test_2 <- "../../../Annotation/GRCh38.p14.primary_assembly.genome.2bit"
  #
  # test_tracks <- inspectorORF::create_tracks(bed_test_2, gtf_test_2, genome_test_2, ms_data_file)
  #
  # test_orf_plot <- inspectorORF::orf_plot(test_tracks, transcript_filter = "ENST00000486256.5", start_position = 1702)
  # test_orf_plot
})

test_that("get ORF nucleotide sequence",
{
  ms_data_file <- "../../../Annotation/control_MS_evidence.bed"

  bed_test_2 <- "../../../Annotation/test2.bed"
  bed_test <- "../../../Annotation/test.bed"

  gtf_test <- "../../../Annotation/gencode.v30.annotation.gtf"
  gtf_test_2 <- "../../../Annotation/gencode.v44.annotation.gtf"

  genome_test <- "../../../Annotation/GRCh38.p12.primary_assembly.genome.2bit"
  genome_test_2 <- "../../../Annotation/GRCh38.p14.primary_assembly.genome.2bit"

  test_tracks <- inspectorORF::create_tracks(bed_test_2, gtf_test_2, genome_test_2, ms_data_file)

  orf_sequence <- inspectorORF::get_orf_nt_seq(test_tracks, transcript_filter = "ENST00000486256.5", start_position = 1702)

  expect_equal(orf_sequence, RNAStringSet(paste0("AUGUUUGCAGAAAUGGAAAUCAUUGGUCAGUUUAACCUGGGAUUUAUAAUAACCAAACUGAAUGAGGAUAUCUUCAUAGUGGACCAGCAUGC",
                                                 "CACGGACGAGAAGUAUAACUUCGAGAUGCUGCAGCAGCACACCGUGCUCCAGGGGCAGAGGCUCAUAGCACCUCAGACUCUCAACUUAACUG",
                                                 "CUGUUAAUGAAGCUGUUCUGAUAGAAAAUCUGGAAAUAUUUAGAAAGAAUGGCUUCGAUUUUGUUAUCGAUGAAAAUGCUCCAGUCACUGAA",
                                                 "AGGGCUAAACUGAUUUCCUUGCCAACUAGUAAAAGCUGGACCUUCGGACCCCAGGACGUCGAUGAACUGAUCUUCAUGCUGAGCGACAGCCC",
                                                 "UGGGGUCAUGUGCCGGCCUUCCCGAGUCAAGCAGAUGUUUGCCUCCAGAGCCUGCCGGAAGUCGGUGAUGAUUGGGACUGCUCUUAACACAA",
                                                 "GCGAGAUGAAGAAACUGAUCACCCACAUGGGGGAGAUGGACCACCCCUGGAACUGUCCCCAUGGAAGGCCAACCAUGAGACACAUCGCCAAC",
                                                 "CUGGGUGUCAUUUCUCAGAACUGA")))
})

test_that("get amino acid ORF sequence",
{
  ms_data_file <- "../../../Annotation/control_MS_evidence.bed"

  bed_test_2 <- "../../../Annotation/test2.bed"
  bed_test <- "../../../Annotation/test.bed"

  gtf_test <- "../../../Annotation/gencode.v30.annotation.gtf"
  gtf_test_2 <- "../../../Annotation/gencode.v44.annotation.gtf"

  genome_test <- "../../../Annotation/GRCh38.p12.primary_assembly.genome.2bit"
  genome_test_2 <- "../../../Annotation/GRCh38.p14.primary_assembly.genome.2bit"

  test_tracks <- inspectorORF::create_tracks(bed_test_2, gtf_test_2, genome_test_2, ms_data_file)

  orf_aa_sequence <- inspectorORF::get_orf_aa_seq(test_tracks, transcript_filter = "ENST00000486256.5", start_position = 1702)

  expect_equal(orf_aa_sequence, AAStringSet(paste0("MFAEMEIIGQFNLGFIITKLNEDIFIVDQHATDEKYNFEMLQQHTVLQGQRLIAP",
                                                    "QTLNLTAVNEAVLIENLEIFRKNGFDFVIDENAPVTERAKLISLPTSKSWTFGPQ",
                                                    "DVDELIFMLSDSPGVMCRPSRVKQMFASRACRKSVMIGTALNTSEMKKLITHMGE",
                                                    "MDHPWNCPHGRPTMRHIANLGVISQN*")))
})
