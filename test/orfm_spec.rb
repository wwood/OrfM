require 'rspec'
require 'bio-commandeer'
require 'tempfile'

orfm = File.expand_path(File.dirname(__FILE__) + '/../orfm')

RSpec.configure do |config|
  config.expect_with :rspec do |c|
    c.syntax = [:should, :expect]
  end
end

describe "orfm" do
  it "should run with defaults on a single file" do
    input = %w(>638202197:1-99 ATGGATGCTGAAAAAAGATTGTTCTTAAAGGCATTAAAGGAAAAGTTTGAAGAAGACCCAAGAGAAAAATACACTAAGTTCTATGTCTTTGGCGGATGG).join("\n")
    expected = %w(>638202197:1-99_1_1_1 MDAEKRLFLKALKEKFEEDPREKYTKFYVFGGW).join("\n")+"\n"
    Bio::Commandeer.run("#{orfm}", :stdin => input).should == expected
  end

  it "should work with a longer -m " do
    input = ['>638202197 NP_247840 methyl coenzyme M reductase I, subunit alpha (mcrA) [Methanocaldococcus jannaschii DSM 2661: NC_000909] (+)strand',
<<EOS
ATGGATGCTGAAAAAAGATTGTTCTTAAAGGCATTAAAGGAAAAGTTTGA
AGAAGACCCAAGAGAAAAATACACTAAGTTCTATGTCTTTGGCGGATGGA
GACAGTCAGCAAGAAAAAGAGAATTCGTTGAGGCAGCACAAAAATTAATT
GAGAAGAGAGGAGGAATTCCATTTTACAACCCAGATATTGGAGTTCCATT
GGGGCAGAGAAAATTAATGCCTTACAAAGTTTCAAATACAGATGCAATTG
TTGAAGGGGATGACTTACACTTCATGAACAACGCTGCAATGCAGCAGTTC
TGGGATGACATAAGAAGAACAGTTATCGTTGGGATGGATACAGCTCACGC
TGTTCTTGAAAAGAGATTGGGGGTAGAGGTTACTCCAGAAACAATTAATG
AATACATGGAAACTATTAACCACGCTCTCCCAGGAGGAGCAGTTGTTCAG
GAGCACATGGTTGAGGTCCACCCAGCATTAGTCTGGGACTGTTACGCTAA
GATATTCACTGGAGATGACGAATTAGCAGATGAGATTGACAAGAGGTTCT
TAATTGACATTAACAAGTTGTTCCCAGAAGAGCAGGCAGAACAAATCAAG
AAGGCAATCGGTAAGAGAACATACCAAGTTTCAAGAGTTCCAACATTAGT
CGGTAGAGTTTGTGATGGGGGAACAATAGCAAGATGGAGTGCTATGCAGA
TTGGAATGTCATTCATTACAGCTTACAAGCTCTGTGCTGGGGAGGCAGCA
ATTGCTGACTTCTCATACGCTGCAAAGCACGCTGATGTCATTCAGATGGC
TTCATTCTTGCCAGCAAGAAGAGCAAGAGGGCCAAATGAACCAGGAGGTA
TCTTCTTCGGAGTCTTGGCAGATATTGTTCAAACATCAAGAGTTTCAGAT
GACCCAGTTGAACAGTCATTAGAGGTTGTTGCTGCTGGGGCTATGTTGTA
TGACCAAATCTGGTTAGGAGGATACATGTCTGGAGGAGTCGGATTTACAC
AGTATGCTACAGCAACCTATACAGATGACATCTTGGATGACTTCTCATAC
TACGGATATGACTACATAACCAAGAAATATGGAGGATGCAACAGCGTAAA
ACCAACAATGGATGTTGTTGAAGATATTGCTACTGAAGTAACTTTATATG
GTTTAGAGCAGTATGACACCTTCCCAGCATTGTTAGAAGACCACTTCGGA
GGTTCCCAAAGAGCAGGGGTTACAGCTGCTGCAGCAGGTATTACAACTGC
ATTAGCTACAGGAAACTCAAACGCTGGAGTTAACGGATGGTATCTAAGCC
AGATATTGCACAAAGAATACCACAGCAGATTAGGATTCTATGGTTATGAC
TTACAAGACCAGTGTGGAGCAGCCAACTCATTATCATTCAGAAACGATGA
AGGTTCCCCATTAGAATTGAGAGGGCCTAACTATCCAAACTACGCAATGA
ACGTTGGTCACCAAGGAGAATATGCTGGAATTACACAGGCTGCACACTCA
GCAAGAGGAGACGCATTTGCATTGAACCCATTAATTAAGGTTGCATTTGC
AGACCCATCATTAGTCTTTGACTTCACACATCCAAGAAAAGAGTTTGCAA
GAGGTGCTTTAAGAGAATTCGAGCCAGCTGGAGAAAGAGATCCAATCATC
CCAGCTCACTAA
EOS
      ].join("\n")
    outname = '>638202197_1_1_1 NP_247840 methyl coenzyme M reductase I, subunit alpha (mcrA) [Methanocaldococcus jannaschii DSM 2661: NC_000909] (+)strand'
    outseq = <<EOS
MDAEKRLFLKALKEKFEEDPREKYTKFYVFGGWRQSARKREFVEAAQKLIEKRGGIPFYN
PDIGVPLGQRKLMPYKVSNTDAIVEGDDLHFMNNAAMQQFWDDIRRTVIVGMDTAHAVLE
KRLGVEVTPETINEYMETINHALPGGAVVQEHMVEVHPALVWDCYAKIFTGDDELADEID
KRFLIDINKLFPEEQAEQIKKAIGKRTYQVSRVPTLVGRVCDGGTIARWSAMQIGMSFIT
AYKLCAGEAAIADFSYAAKHADVIQMASFLPARRARGPNEPGGIFFGVLADIVQTSRVSD
DPVEQSLEVVAAGAMLYDQIWLGGYMSGGVGFTQYATATYTDDILDDFSYYGYDYITKKY
GGCNSVKPTMDVVEDIATEVTLYGLEQYDTFPALLEDHFGGSQRAGVTAAAAGITTALAT
GNSNAGVNGWYLSQILHKEYHSRLGFYGYDLQDQCGAANSLSFRNDEGSPLELRGPNYPN
YAMNVGHQGEYAGITQAAHSARGDAFALNPLIKVAFADPSLVFDFTHPRKEFARGALREF
EPAGERDPIIPAH
EOS
    outseq.gsub!("\n",'')

    Bio::Commandeer.run("#{orfm} -m 300", :stdin => input).should == "#{outname}\n#{outseq}\n"
  end

  it 'should toy example with internal frame 2 ORF' do
    input = %w(>eg AATGTGAA).join("\n")
    expected = <<EOS
>eg_2_2_1
M
>eg_1_1_2
NV
>eg_3_3_3
CE
>eg_1_4_4
HI
>eg_2_5_5
SH
>eg_3_6_6
FT
EOS
    Bio::Commandeer.run("#{orfm} -m3", :stdin => input).should == expected
  end

  it 'should be able to handle n characters' do
    input = %w(>eg TTAANA).join("\n")
    expected = %w(>eg_1_1_1 LX).join("\n")+"\n";
    Bio::Commandeer.run("#{orfm} -m6", :stdin => input).should == expected
  end

  it 'should be able to handle lower case characters' do
    input = %w(>eg TTAAaA).join("\n")
    expected = %w(>eg_1_1_1 LK).join("\n")+"\n";
    Bio::Commandeer.run("#{orfm} -m6", :stdin => input).should == expected
  end

  it 'should exit with non-zero status when -l < -m' do
    input = %w(>eg TTAAaA).join("\n")
    expect {
      Bio::Commandeer.run("#{orfm} -l 3 -m6", :stdin => input).should == expected
      }.to raise_exception
  end

  it 'should stop when it runs out of -l' do
    input = %w(>eg TTAANAGGGGGGGGGG).join("\n")
    expected = %w(>eg_1_1_1 LX).join("\n")+"\n";
    Bio::Commandeer.run("#{orfm} -m6 -l6", :stdin => input).should == expected

    expected = %w(>eg_1_1_1 LX >eg_2_5_2 XL).join("\n")+"\n";
    Bio::Commandeer.run("#{orfm} -m6 -l7", :stdin => input).should == expected

    expected = %w(>eg_1_1_1 LX >eg_3_3_2 XR >eg_2_5_3 XL >eg_3_6_4 PX).join("\n")+"\n";
    Bio::Commandeer.run("#{orfm} -m6 -l8", :stdin => input).should == expected

    expected = %w(>eg_1_1_1
       LXGGG
       >eg_5_2_2
       XGGG
       >eg_3_3_3
       XRGG
       >eg_4_4_4
       PPPX
       >eg_2_5_5
       PPPXL
       >eg_3_6_6
       PPPX).join("\n")+"\n";
    Bio::Commandeer.run("#{orfm} -m6 -l16", :stdin => input).should == expected
    Bio::Commandeer.run("#{orfm} -m6 -l80", :stdin => input).should == expected
  end

  it 'should require versions as you would expect' do
    vers = Bio::Commandeer.run("#{orfm} -v").strip.split(' ')[2].split('.').collect{|v| v.to_i}

    input = %w(>eg TTAAaA).join("\n")
    expect {Bio::Commandeer.run("#{orfm} -r #{vers[0]+1}.#{vers[1] }.#{vers[2] }", :stdin => input).should == expected}.to raise_exception
    expect {Bio::Commandeer.run("#{orfm} -r #{vers[0] }.#{vers[1]+1}.#{vers[2] }", :stdin => input).should == expected}.to raise_exception
    expect {Bio::Commandeer.run("#{orfm} -r #{vers[0] }.#{vers[1] }.#{vers[2]+1}", :stdin => input).should == expected}.to raise_exception
    Bio::Commandeer.run("#{orfm} -r #{vers[0] }.#{vers[1] }.#{vers[2] }", :stdin => input).should == ''
    Bio::Commandeer.run("#{orfm} -r #{vers[0] }.#{vers[1]-1}.#{vers[2] }", :stdin => input).should == ''
    Bio::Commandeer.run("#{orfm} -r #{vers[0] }.#{vers[1]-1}.#{vers[2]+1}", :stdin => input).should == ''
  end

  it 'should accept bad codons' do
    input = %w(>eg GYATCATAGGCCAGCCGCTGTCCAGATGCACCGGTTCATCTGCGTCAGACGACGATCTTCACCCGGTAACCCCCGCCGATCACCAGATACTCGGCCTCCCCGCGCAAAGGATGGCTCGGCAGCAGATCGTTGAAGAACAGGAGCTTCACGACCGGAACTTCGGTTTCGAGGATATAGGCACCGAACTGCCCCGCCTGCGTCCGGTCAGCGGAGAAAGAAACGATGTTGTTGAGACGCACGAGGATTTCCCGTCCGTTGCCGGCCCC).join("\n")
    expected = %w(>eg_1_1_1
       XS
       >eg_2_2_2
       XH
       >eg_3_3_3
       II
       >eg_2_5_4
       MX
       >eg_3_6_5
       YD).join("\n")+"\n";
    Bio::Commandeer.run("#{orfm} -m6 -l8", :stdin => input).should == expected
  end

  it 'should give the right transcripts' do
    input = %w(>eg AATGTGAA).join("\n")
    expected = %w(>eg_2_2_1
        ATG
        >eg_1_1_2
        AATGTG
        >eg_3_3_3
        TGTGAA
        >eg_1_4_4
        CACATT
        >eg_2_5_5
        TCACAT
        >eg_3_6_6
        TTCACA)
    Tempfile.open("orfm_testing") do |t|
      Bio::Commandeer.run("#{orfm} -m3 -t #{t.path}", :stdin => input)
      File.open(t.path).read.split("\n").should == expected
    end
  end

  it 'should handle odd chars in transcript output' do
    input = %w(>eg GYATCATAGGCCAGCCGCTGTCCAGATGCACCGGTTCATCTGCGTCAGACGACGATCTTCACCCGGTAACCCCCGCCGATCACCAGATACTCGGCCTCCCCGCGCAAAGGATGGCTCGGCAGCAGATCGTTGAAGAACAGGAGCTTCACGACCGGAACTTCGGTTTCGAGGATATAGGCACCGAACTGCCCCGCCTGCGTCCGGTCAGCGGAGAAAGAAACGATGTTGTTGAGACGCACGAGGATTTCCCGTCCGTTGCCGGCCCC).join("\n")
    expected = %w(>eg_1_1_1
       GYATCA
       >eg_2_2_2
       YATCAT
       >eg_3_3_3
       ATCATA
       >eg_2_5_4
       ATGATN
       >eg_3_6_5
       TATGAT)
    Tempfile.open("orfm_testing") do |t|
      Bio::Commandeer.run("#{orfm} -t #{t.path} -m6 -l8", :stdin => input)
      File.open(t.path).read.split("\n").should == expected
    end
  end

  it 'should accept alternative codon tables' do
    input = %w(>eg AATGTGAA).join("\n")
    expected = %w(
        >eg_1_1_1
        NV
        >eg_2_2_2
        MW
        >eg_3_3_3
        CE
        >eg_1_4_4
        HI
        >eg_2_5_5
        SH
        >eg_3_6_6
        FT)
    Bio::Commandeer.run("#{orfm} -m3 -c4", :stdin => input).split("\n").should == expected
  end
end
