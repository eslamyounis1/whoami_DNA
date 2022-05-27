import streamlit as st
import neatbio.sequtils as utils
import io
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import random
from Bio.Seq import Seq
from Bio import Seq
from Bio import SeqIO
from collections import Counter
from PIL import Image
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from bs4 import BeautifulSoup
import requests
import csv
from itertools import zip_longest
matplotlib.use("Agg")


def delta(x, y):
    return 0 if x == y else 1


def M(seq1, seq2, i, j, k):
    return sum(delta(x, y) for x, y in zip(seq1[i:i + k], seq2[j:j + k]))


def makeMatrix(seq1, seq2, k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1, seq2, i, j, k) for j in range(m - k + 1)] for i in range(n - k + 1)]


def plotMatrix(M, t, seq1, seq2, nonblank=chr(0x25A0), blank=' '):
    print(' |' + seq2)
    print('-' * (2 + len(seq2)))
    for label, row in zip(seq1, M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)


def dotplot(seq1, seq2, k=1, t=1):
    M = makeMatrix(seq1, seq2, k)
    plotMatrix(M, t, seq1, seq2)  # experiment with character choice


# Convert to Fxn
def dotplotx(seq1, seq2):
    plt.imshow(np.array(makeMatrix(seq1, seq2, 1)))
    # on x-axis list all sequences of seq 2
    xt = plt.xticks(np.arange(len(list(seq2))), list(seq2))
    # on y-axis list all sequences of seq 1
    yt = plt.yticks(np.arange(len(list(seq1))), list(seq1))
    plt.show()


def gc_content(seq):
    result = float(str(seq).count('G') + str(seq).count('C')) / len(seq) * 100
    return result


def at_content(seq):
    result = float(str(seq).count('A') + str(seq).count('T')) / len(seq) * 100
    return result


def insertion_mutations(sequence):
    bases = 'ATCG'
    insert_bases = ''
    seq = ''
    range_mutation = random.randint(50, 100)

    for rand in range(range_mutation):
        random_number = random.randint(0, 1)
        if random_number == 0:
            for i in range(100):
                number = random.randint(0, 3)
                insert_bases += bases[number]

            index = random.randint(0, len(sequence) - 1)
            seq = sequence[:index] + insert_bases + sequence[index:]
        else:
            number = random.randint(0, 3)
            insert_bases += bases[number]
            index = random.randint(0, len(sequence) - 1)
            seq = sequence[:index] + insert_bases + sequence[index:]

    with st.expander("See Sequence"):
        st.write(seq)

    st.download_button("Download Insertion File", str(seq), file_name="insertion_mutations.txt")
    return seq


def deletion_mutations(sequence):
    seq = Seq.MutableSeq(sequence)

    range_mutation = random.randint(1000, 4000)
    for rand in range(range_mutation):
        number = random.randint(0, len(seq) - 1)
        seq.remove(seq[number])

    with st.expander("See Sequence"):
        st.write(seq)

    st.download_button("Download Deletion File", str(seq), file_name="deletion_mutations.txt")
    return seq


def substitution_mutations(sequence):
    seq = Seq.MutableSeq(sequence)
    bases = 'ATCG'
    number = random.randint(100, 300)
    for i in range(number):
        number_bases = random.randint(0, 3)
        number_seq = random.randint(0, len(seq))
        seq[number_seq] = bases[number_bases]

    with st.expander("See Sequence"):
        st.write(seq)

    st.download_button("Download Substitution File", str(seq), file_name="substitution_mutations.txt")
    return seq


def visualize_together(normal, mutation):
    sorted(normal)
    sorted(mutation)
    col1, col2, col3 = st.columns(3)

    with col1:
        st.header("violinplot")
        plt.violinplot(list(normal.values()))
        plt.violinplot(list(mutation.values()))
        # plt.legend(['normal', 'mutation'])
        st.set_option('deprecation.showPyplotGlobalUse', False)
        st.pyplot()

    with col2:
        st.header("plot")

        plt.plot(list(normal.values()))
        plt.plot(list(mutation.values()), linestyle='dashdot')
        plt.legend(['normal', 'mutation'])
        st.set_option('deprecation.showPyplotGlobalUse', False)
        st.pyplot()

    with col3:
        st.header("boxplot")

        plt.boxplot(list(normal.values()))
        plt.boxplot(list(mutation.values()))
        # plt.legend(['normal', 'mutation'])
        st.set_option('deprecation.showPyplotGlobalUse', False)
        st.pyplot()


def count_dict(Translate):
    dict = {}
    for base in Translate.replace('*', ''):
        if base in dict:
            dict[base] += 1
        else:
            dict[base] = 0

    return dict


def main():
    # st.snow()
    image = Image.open('./logo/logo.png')

    st.image(image, use_column_width=True)

    # """A Simple Bioinformatics App"""
    st.title("WHOAMI - Sequencing App")

    menu = ["Intro", "DNA Sequence", "DotPlot", "Alignment", "DNA Mutation", "Blast", "About"]
    choice = st.sidebar.selectbox("Select Option", menu)

    if choice == "Intro":
        st.subheader("Python Project -Bioinformatics")

    elif choice == "Alignment":

        seq1 = st.text_input("Enter the 1st Sequence").upper()
        seq2 = st.text_input("Enter the 2nd Sequence").upper()
        # # Finding similarities
        alignments = pairwise2.align.globalxx(seq1, seq2)

        # Showing results

        if st.button("Align"):

            with st.container():
                with st.expander("See Sequence"):
                    for alignment in alignments:
                        st.text(format_alignment(*alignment))


    elif choice == "DNA Mutation":

        st.subheader("DNA Mutation")
        seq_file_mutation = st.file_uploader("Upload FASTA File", type=["fasta", "fa", "fastq"])

        if seq_file_mutation is not None:
            byte_str3 = seq_file_mutation.read()
            text_obj3 = byte_str3.decode('UTF-8')
            dna_record = SeqIO.read(io.StringIO(text_obj3), "fasta")
            # st.write(pep_seq)
            dna_seq = dna_record.seq

            details = st.radio("Details", ("Description", "Sequence"))
            if details == "Description":
                st.write(dna_record.description)
            elif details == "Sequence":
                with st.expander("See Sequence"):
                    st.write(dna_record.seq)
                st.write("--------------------------------------------------------")
                with st.container():
                    if st.checkbox("Insertion Mutation", value=True):
                        mutate_seq = insertion_mutations(dna_record.seq)
                        if st.button("Insertion Plot"):
                            st.write('DNA Sequence:')
                            dna_mute_freq = count_dict(mutate_seq)
                            dna_norm_freq = count_dict(dna_record.seq)
                            visualize_together(dna_norm_freq, dna_mute_freq)

                            st.write('amino acids:')
                            dna_mute_freq = count_dict(mutate_seq.translate())
                            dna_norm_freq = count_dict(dna_record.seq.translate())
                            visualize_together(dna_norm_freq, dna_mute_freq)
                st.write("--------------------------------------------------------")
                with st.container():
                    if st.checkbox("Substitution Mutation", value=True):
                        mutate_seq = substitution_mutations(dna_record.seq)
                        if st.button("Substitution Plot"):
                            st.write('DNA Sequence:')
                            dna_mute_freq = count_dict(mutate_seq)
                            dna_norm_freq = count_dict(dna_record.seq)
                            visualize_together(dna_norm_freq, dna_mute_freq)

                            st.write('amino acids:')
                            dna_mute_freq = count_dict(mutate_seq.translate())
                            dna_norm_freq = count_dict(dna_record.seq.translate())
                            visualize_together(dna_norm_freq, dna_mute_freq)

                st.write("--------------------------------------------------------")
                with st.container():
                    if st.checkbox("Deletion Mutation", value=True):
                        mutate_seq = deletion_mutations(dna_record.seq)
                        if st.button("Deletion Plot"):
                            st.write('DNA Sequence:')
                            dna_mute_freq = count_dict(mutate_seq)
                            dna_norm_freq = count_dict(dna_record.seq)
                            visualize_together(dna_norm_freq, dna_mute_freq)

                            st.write('amino acids:')
                            dna_mute_freq = count_dict(mutate_seq.translate())
                            dna_norm_freq = count_dict(dna_record.seq.translate())
                            visualize_together(dna_norm_freq, dna_mute_freq)






    elif choice == "DNA Sequence":
        st.subheader("DNA Sequence Analysis")

        seq_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa", "fastq"])

        if seq_file is not None:
            byte_str = seq_file.read()
            text_obj = byte_str.decode('UTF-8')
            dna_record = SeqIO.read(io.StringIO(text_obj), "fasta")
            # st.write(pep_seq)
            dna_seq = dna_record.seq

            details = st.radio("Details", ("Description", "Sequence"))
            if details == "Description":
                st.write(dna_record.description)
            elif details == "Sequence":
                with st.expander("See Sequence"):
                    st.write(dna_record.seq)

                # Nucleotide Frequencies
                st.subheader("Nucleotide Frequency")
                dna_freq = Counter(dna_seq)
                with st.expander("See Sequence"):
                    st.write(dna_freq)

                guanine_color = st.color_picker("Guanine Color")
                thymine_color = st.color_picker("ThymineColor")
                adenine_color = st.color_picker("Adenine Color")
                cytosine_color = st.color_picker("Cytosine Color")

                if st.button("Nucleotide Plot"):
                    barlist = plt.bar(dna_freq.keys(), dna_freq.values())
                    barlist[0].set_color(guanine_color)
                    barlist[1].set_color(thymine_color)
                    barlist[2].set_color(adenine_color)
                    barlist[3].set_color(cytosine_color)
                    st.set_option('deprecation.showPyplotGlobalUse', False)
                    st.pyplot()
                st.subheader("DNA Composition Scores")
                gc_score = utils.gc_content(str(dna_seq))
                at_score = utils.at_content(str(dna_seq))
                st.json({"GC Content": gc_score, "AT Content": at_score})

                # Nucleotide count search
                nt_count = st.text_input("Enter Nucleotide Here: ", "Type Nucleotide Alphabet")
                st.write(("Number of {} Nucleotide is ::{}".format(nt_count, str(dna_seq).count(nt_count))))

                # Protein Synthesis - Transcription - Translation
                st.subheader("Protein Synthesis")
                p1 = dna_seq.translate()
                aa_freq = Counter(str(p1))

                if st.checkbox("Transcription"):
                    transcription = dna_seq.transcribe()
                    with st.expander("See Sequence"):
                        st.write(transcription)
                    st.download_button("Download Transcription", str(transcription),
                                       file_name='seq_transcription.txt')
                if st.checkbox("Translation"):
                    translation = dna_seq.translate()
                    with st.container():
                        with st.expander("See Sequence"):
                            st.write(translation)
                        st.download_button("Download Transcription", str(translation),
                                           file_name='seq_translation.txt')
                if st.checkbox("Complement"):
                    complement = dna_seq.complement()
                    with st.container():
                        with st.expander("See Sequence"):
                            st.write(complement)
                        st.download_button("Download Transcription", str(complement),
                                           file_name='seq_complement.txt')
                if st.checkbox("AA Frequency"):
                    with st.container():
                        with st.expander("See Sequence"):
                            st.write(aa_freq)

                if st.checkbox("AA Plot"):
                    aa_color = st.color_picker("Pick An Amino Acid Color")
                    # barlist = plt.bar(aa_freq.keys(), aa_freq.values())
                    # barlist[0].set_color(aa_color)
                    plt.bar(aa_freq.keys(), aa_freq.values(), color=aa_color)
                    st.set_option('deprecation.showPyplotGlobalUse', False)
                    st.pyplot()

                if st.checkbox("Full Amino Acid Name"):
                    aa_name = str(p1).replace("*", "")
                    aa3 = utils.convert_1to3(aa_name)
                    with st.container():
                        with st.expander("See Sequence"):
                            with st.container():
                                st.write(aa_name)
                        st.write("--------------------------------------------------------")
                        with st.expander("See Sequence"):
                            with st.container():
                                st.write(aa3)

                        st.write("--------------------------------------------------------")
                        get_a = utils.get_acid_name(aa3)
                        with st.expander("See Sequence"):
                            st.write(get_a)

                        st.download_button("Download Full Amino Acid Name", str(get_a), file_name="AA Full Name.txt")

                # Top Most Common Amino
    elif choice == "DotPlot":

        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file1 = st.file_uploader("Upload 1st FASTA File", type=["fasta", "fa", "fastq"])
        seq_file2 = st.file_uploader("Upload 2nd FASTA File", type=["fasta", "fa", "fastq"])

        if seq_file1 and seq_file2 is not None:
            byte_str2 = seq_file1.read()
            text_obj1 = byte_str2.decode('UTF-8')
            dna_record1 = SeqIO.read(io.StringIO(text_obj1), "fasta")
            dna_seq1 = dna_record1.seq

            byte_str3 = seq_file2.read()
            text_obj2 = byte_str3.decode('UTF-8')
            dna_record2 = SeqIO.read(io.StringIO(text_obj2), "fasta")
            dna_seq2 = dna_record2.seq

            details = st.radio("Details", ("Description", "Sequence"))
            if details == "Description":
                st.write(dna_record1.description)
                st.write("--------------------------------------------------------")
                st.write(dna_record2.description)
            elif details == "Sequence":

                with st.expander("See Sequence of sample 1"):
                    st.write(dna_record1.seq)
                st.write("--------------------------------------------------------")

                with st.expander("See Sequence of sample 2"):
                    st.write(dna_record2.seq)

            cus_limit = st.number_input("Select Max number of Nucleotide", 10, 200, 50)
            if st.button("Dot Plot"):
                st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                dotplotx(dna_seq1[0:cus_limit], dna_seq2[0:cus_limit])
                st.set_option('deprecation.showPyplotGlobalUse', False)

                st.pyplot()

    elif choice == "Blast":
        st.subheader("PSI-BLAST")
        # prot_input = st.text_input("Enter protein sequence")
        html_code = '<a href="https://www.ebi.ac.uk/Tools/sss/psiblast/" target="_blank">GO TO NCBI</a>'
        st.markdown(html_code, unsafe_allow_html=True)
        url_text = st.text_input("Enter URL")

        if st.button("Protein Result"):
            results = requests.get(url_text)

            ID = []
            length = []
            score = []
            identity = []
            positives = []
            E = []

            src = results.content
            soup = BeautifulSoup(src, 'lxml')

            for i in range(2):
                Alignment_results = soup.find_all('tr', {'id': 'alignment_' + str(i)})
                Alignment = Alignment_results[0].text.split('\n')

                ID.append(Alignment[6])
                length.append(Alignment[12])
                score.append(Alignment[13])
                identity.append(Alignment[14])
                positives.append(Alignment[15])
                E.append(Alignment[16])

            file_lists = [ID, length, score, identity, positives]
            export = zip_longest(*file_lists)

            alignment = open('.\\alignment.csv', 'w+')
            wr = csv.writer(alignment)
            wr.writerow(['ID', 'length', 'score', 'identity', 'positives', 'E'])
            wr.writerows(export)

            st.success('Created Successfully')
            st.download_button("Download File", str(alignment), file_name="alignment.csv")


        if st.button("Visualize Similarity"):
            with open('.\\alignment.csv', 'r') as alignment:
                file = alignment.read()
                file = file.split('\n')
                lines = []
                for i in range(len(file)):
                    if len(file[i]) > 0:
                        lines.append(file[i].split(','))

                score = []
                identity = []
                positive = []
                IDs = []
                for i in range(1, len(lines)):
                    IDs.append(lines[i][0])
                    score.append(lines[i][2])
                    identity.append(lines[i][3])
                    positive.append(lines[i][4])

                st.set_option('deprecation.showPyplotGlobalUse', False)
                with st.container():
                    col1, col2 = st.columns(2)
                    with col1:
                        ex = []
                        plt.title('score')
                        for i in range(len(score)):
                            ex.append(0.1)
                        plt.pie(score, labels=IDs, autopct="%2.1f%%", explode=ex)
                        st.pyplot()
                    with col2:
                        ex = []
                        plt.title('identity')
                        for i in range(len(identity)):
                            ex.append(0.1)
                        plt.pie(identity, labels=IDs, autopct="%2.1f%%", explode=ex)
                        st.pyplot()


    elif choice == "About":

        st.subheader("About")
        st.text('''
                        Team Members:
                               1- Eslam Younis
                               2- Maria Sameh
                               3- Ester elzek
                               4- Mira Rafat
                               5- Merihan Yousef
                               6- Steven Waheed''')


if __name__ == '__main__':
    main()

