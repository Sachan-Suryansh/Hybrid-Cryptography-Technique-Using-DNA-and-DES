import os
import sys

from time import time
import ast
import utils
from utils import *
import string
from time import time

from PyQt5 import QtCore
from PyQt5.QtCore import QCoreApplication, QTimer, pyqtSlot
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QMainWindow
from PyQt5.uic import loadUi

from des import DesKey


class DNA_Encryption_Functions():
    """
    DNA Genetic Encryption
    """

    # number of rounds the algorithm is run, chosen randomly
    rounds_no = None

    chromosome_length = None

    # the key used in the decryption process
    decryption_key = None

    def set_globals(self):
        global rounds_no
        global decryption_key
        # it is better to be odd random
        rounds_no = random.randrange(3, 12, 2)
        decryption_key = ""

    def encrypt_key(self, data, key):
        """
        Encrypt data with key: data XOR key.
        """

        # repeat key ONLY if data is longer than key and encrypt
        if len(data) > len(key):
            factor = int(len(data) / len(key))
            key += key * factor

            return bitxor(data, key)

        return bitxor(data, key)

    def reshape(self, dna_sequence):
        """
        Generate chromosome population.
        :param dna_sequence: a string sequence of DNA bases
        :return: an array of chromosomes, chromosome population
        """
        global chromosome_length
        global decryption_key

        # choose population size and chromosome length
        divs = divisors(len(dna_sequence))
        chromosome_no = divs[random.randint(0, len(divs) - 1)]
        chromosome_length = int(len(dna_sequence) / chromosome_no)
        chromosomes = []

        decryption_key += reshape_del + str(chromosome_length) + reshape_del

        # retrieve the population
        for i in range(0, len(dna_sequence), chromosome_length):
            chromosomes.append(dna_sequence[i:i + chromosome_length])

        return chromosomes

    def reverse_reshape(self, population):
        # convert the chromosome population back to DNA sequence
        return "".join(population)

    def rotate_crossover(self, population):
        """
        Rotate every chromosome in population left / right according to probability p.
        """
        global chromosome_length
        global decryption_key

        new_population = []

        decryption_key += rotate_crossover_del

        # predefined rotation value, varied every round
        rotation_offset = random.randint(1, chromosome_length)

        decryption_key += rotation_offset_del + str(rotation_offset) + rotation_offset_del

        decryption_key += rotation_types_del

        for chromosome in population:

            p = random.uniform(0, 1)

            if p > 0.5:
                decryption_key += "right|"
                right_first = chromosome[0: len(chromosome) - rotation_offset]
                right_second = chromosome[len(chromosome) - rotation_offset:]
                new_population.append(right_second + right_first)
            else:
                decryption_key += "left|"
                left_first = chromosome[0: rotation_offset]
                left_second = chromosome[rotation_offset:]
                new_population.append(left_second + left_first)

        decryption_key += rotation_types_del

        decryption_key += rotate_crossover_del

        return new_population

    def single_point_crossover(self, population):
        """
        Combine each two chromosomes in population by using single point crossover.
        """
        global decryption_key

        decryption_key += single_point_crossover_del

        new_population = []
        for i in range(0, len(population) - 1, 2):
            candidate1 = population[i]
            candidate2 = population[i + 1]

            # chromosomes have the same length
            # choose a random point
            length = len(candidate1)
            crossover_point = random.randint(0, length - 1)

            decryption_key += str(crossover_point) + "|"

            offspring1 = candidate2[0: crossover_point] + candidate1[crossover_point:]
            offspring2 = candidate1[0: crossover_point] + candidate2[crossover_point:]
            new_population.append(offspring1)
            new_population.append(offspring2)

        # append last chromosome if odd population size
        if len(population) % 2 == 1:
            new_population.append(population[len(population) - 1])

        decryption_key += single_point_crossover_del

        return new_population

    def crossover(self, population):
        global decryption_key

        # choose crossover type according to p
        p = random.uniform(0, 1)

        if p < 0.33:
            decryption_key += crossover_type_del + "rotate_crossover" + crossover_type_del
            return self.rotate_crossover(population)
        elif p >= 0.33 and p < 0.66:
            decryption_key += crossover_type_del + "single_point_crossover" + crossover_type_del
            return self.single_point_crossover(population)
        else:
            decryption_key += crossover_type_del + "both" + crossover_type_del
            population = self.rotate_crossover(population)
            return self.single_point_crossover(population)

    def complement(self, chromosome, point1, point2):
        """
        Flip chromosome bits between point1 and point2.
        """
        new_chromosome = ""

        for i in range(len(chromosome)):
            if i >= point1 and i <= point2:
                if chromosome[i] == '0':
                    new_chromosome += '1'
                else:
                    new_chromosome += '0'
            else:
                new_chromosome += chromosome[i]

        return new_chromosome

    def alter_dna_bases(self, bases):
        """
        Alter DNA bases to another one randomly.(e.g. C->G and A->T and viceversa)
        """
        alter_dna_table = {}

        for _ in range(2):
            # choose one randomly then remove it from list
            base1 = bases[random.randint(0, len(bases) - 1)]
            bases.remove(base1)

            # choose one randomly then remove it from list
            base2 = bases[random.randint(0, len(bases) - 1)]
            bases.remove(base2)

            # assign the first to the other
            alter_dna_table[base1] = base2
            alter_dna_table[base2] = base1

        return alter_dna_table

    def mutation(self, population):
        """
        Apply mutation operator by using "complement" and "alter_dna_bases"
        """
        global decryption_key

        bases = ['A', 'C', 'G', 'T']
        alter_dna_table = self.alter_dna_bases(bases)

        decryption_key += mutation_table_del + str(alter_dna_table) + mutation_table_del

        new_population = []
        for chromosome in population:
            decryption_key += chromosome_del

            # apply the complement
            b_chromosome = dna_to_bits(chromosome, utils.dna_base_to_two_bits_table)
            decryption_key += complement_mutation_del
            point1 = random.randint(0, len(b_chromosome) - 1)
            point2 = random.randint(point1, len(b_chromosome) - 1)
            decryption_key += "(%s, %s)" % (point1, point2)
            decryption_key += complement_mutation_del
            b_chromosome = self.complement(b_chromosome, point1, point2)

            # convert each 4 bits in chromosome to two dna bases using four_bits_to_two_dna_base_table
            four_bits_vector = group_bits(b_chromosome, 4)

            last_dna_base = None
            # if the last element is of length 2, don't convert it
            if len(four_bits_vector[len(four_bits_vector) - 1]) == 2:
                last_dna_base = utils.two_bits_to_dna_base_table[four_bits_vector[len(four_bits_vector) - 1]]

                # convert only the 4 bits elements
                four_bits_vector = four_bits_vector[:-1]

            dna_seq = bits_to_dna(four_bits_vector, utils.four_bits_to_two_dna_base_table)
            if last_dna_base is not None:
                dna_seq += last_dna_base

            # and then alter the dna bases between point1 and point2
            decryption_key += alter_mutation_del
            point1 = random.randint(0, len(dna_seq) - 1)
            point2 = random.randint(point1, len(dna_seq) - 1)
            decryption_key += "(%s, %s)" % (point1, point2)
            decryption_key += alter_mutation_del
            new_chromosome = ""
            for i in range(len(dna_seq)):
                if i >= point1 and i <= point2:
                    new_chromosome += alter_dna_table[dna_seq[i]]
                else:
                    new_chromosome += dna_seq[i]

            new_population.append(new_chromosome)

            decryption_key += chromosome_del

        return new_population

    def dna_get(self, text, key):
        global rounds_no
        global decryption_key

        print("\nDNA-GET is running...\n")

        # binarize data and convert it to dna sequence
        b_data1 = binarized_data(text)
        dna_seq = bits_to_dna(b_data1, utils.two_bits_to_dna_base_table)
        # print(dna_seq)

        # there is no need for first reshape like in the pseudocode because my reverse_reshape can work on dna_seq, too
        # i.e. ("".join("ACGT") -> "ACGT")

        b_data2 = dna_seq
        print("Initial DNA sequence:", dna_seq)

        decryption_key += no_rounds_del + str(rounds_no) + no_rounds_del

        # run the algorithm "rounds_no" times
        while rounds_no > 0:
            decryption_key += round_del

            # encrypt data with key after reshaping it back to binary sequence and then convert it back to dna sequence
            b_data2 = bits_to_dna(
                group_bits(
                    self.encrypt_key(dna_to_bits(self.reverse_reshape(b_data2), utils.dna_base_to_two_bits_table),
                                     key)),
                utils.two_bits_to_dna_base_table)
            # print("Encrypted data:", b_data2)

            # create the chromosome population
            b_data2 = self.reshape(b_data2)
            # print("Population data:", b_data2)

            # apply crossover on population
            decryption_key += crossover_del
            b_data2 = self.crossover(b_data2)
            decryption_key += crossover_del
            # print("Population data:", b_data2)

            # apply mutation on population
            decryption_key += mutation_del
            b_data2 = self.mutation(b_data2)
            decryption_key += mutation_del
            # print("Population data:", b_data2)

            rounds_no -= 1

            decryption_key += round_del

        return self.reverse_reshape(b_data2)


class DNA_Decryption_Functions():
    """
    DNA Genetic Decryption
    """

    def encrypt_key(self, data, key):
        """
        Encrypt data with key: data XOR key.
        """

        # repeat key ONLY if data is longer than key and encrypt
        if len(data) > len(key):
            factor = int(len(data) / len(key))
            key += key * factor

            return bitxor(data, key)

        return bitxor(data, key)

    def reshape(self, dna_sequence, reshape_info):
        """
        Generate chromosome population.
        :param dna_sequence: a string sequence of DNA bases
        :param reshape_info: the length of each chromosome
        :return: an array of chromosomes, chromosome population
        """

        chromosome_length = int(reshape_info[0])
        chromosomes = []

        # retrieve the population
        for i in range(0, len(dna_sequence), chromosome_length):
            chromosomes.append(dna_sequence[i:i + chromosome_length])

        return chromosomes

    def reverse_reshape(self, population):
        # convert the chromosome population back to DNA sequence
        return "".join(population)

    def rotate_crossover(self, population, rotate_info):
        """
        Rotate every chromosome in population left / right according to probability p.
        """

        new_population = []

        # get the rotation value
        rotation_offset = int(get_pattern(rotation_offset_del, rotate_info)[0])

        rotations = get_pattern(rotation_types_del, rotate_info)[0].split("|")[:-1]

        for i in range(len(population)):
            chromosome = population[i]

            direction = rotations[i]

            if direction == "left":
                right_first = chromosome[0: len(chromosome) - rotation_offset]
                right_second = chromosome[len(chromosome) - rotation_offset:]
                new_population.append(right_second + right_first)
            elif direction == "right":
                left_first = chromosome[0: rotation_offset]
                left_second = chromosome[rotation_offset:]
                new_population.append(left_second + left_first)

        return new_population

    def single_point_crossover(self, population, single_point_info):
        """
        Combine each two chromosomes in population by using single point crossover.
        """
        crossover_points = [int(p) for p in single_point_info.split("|") if p != '']

        new_population = []
        for i in range(0, len(population) - 1, 2):
            candidate1 = population[i]
            candidate2 = population[i + 1]

            # get the crossover_point
            crossover_point = crossover_points[int(i / 2)]

            offspring1 = candidate2[0: crossover_point] + candidate1[crossover_point:]
            offspring2 = candidate1[0: crossover_point] + candidate2[crossover_point:]
            new_population.append(offspring1)
            new_population.append(offspring2)

        # append last chromosome if odd population size
        if len(population) % 2 == 1:
            new_population.append(population[len(population) - 1])

        return new_population

    def crossover(self, population, crossover_info):
        # get the crossover type
        crossover_type = get_pattern(crossover_type_del, crossover_info)[0]

        if crossover_type == "rotate_crossover":
            rotate_info = get_pattern(rotate_crossover_del, crossover_info)[0]
            return self.rotate_crossover(population, rotate_info)
        elif crossover_type == "single_point_crossover":
            single_point_info = get_pattern(single_point_crossover_del, crossover_info)[0]
            return self.single_point_crossover(population, single_point_info)
        elif crossover_type == "both":
            rotate_info = get_pattern(rotate_crossover_del, crossover_info)[0]
            single_point_info = get_pattern(single_point_crossover_del, crossover_info)[0]
            population = self.single_point_crossover(population, single_point_info)
            return self.rotate_crossover(population, rotate_info)

    def complement(self, chromosome, point1, point2):
        """
        Flip chromosome bits between point1 and point2.
        """
        new_chromosome = ""

        for i in range(len(chromosome)):
            if i >= point1 and i <= point2:
                if chromosome[i] == '0':
                    new_chromosome += '1'
                else:
                    new_chromosome += '0'
            else:
                new_chromosome += chromosome[i]

        return new_chromosome

    def mutation(self, population, mutation_info):
        """
        Apply mutation operator by using "complement" and "alter_dna_bases"
        """

        # extract the alteration table
        alter_dna_table = ast.literal_eval(get_pattern(mutation_table_del, mutation_info[0])[0])

        chromosomes_info = get_pattern(chromosome_del, mutation_info[0])

        new_population = []
        for i in range(len(population)):
            chromosome = population[i]
            chromosome_info = chromosomes_info[i]

            # alter back the dna bases between point1 and point2
            alter_info = get_pattern(alter_mutation_del, chromosome_info)[0]
            point1, point2 = ast.literal_eval(alter_info)
            new_chromosome = ""
            for i in range(len(chromosome)):
                if i >= point1 and i <= point2:
                    new_chromosome += alter_dna_table[chromosome[i]]
                else:
                    new_chromosome += chromosome[i]

            two_bases_vector = group_bases(new_chromosome)

            # last base was not converted using four_bits_to_two_dna_base_table
            # convert it to bits using dna_base_to_two_bits_table
            last_two_bits = None
            if len(new_chromosome) % 2 == 1:
                last_two_bits = utils.dna_base_to_two_bits_table[new_chromosome[-1]]

                two_bases_vector = two_bases_vector[:-1]

            bits_seq = dna_to_bits(two_bases_vector, utils.two_dna_base_to_four_bits_table)

            if last_two_bits is not None:
                bits_seq += last_two_bits

            complement_info = get_pattern(complement_mutation_del, chromosome_info)[0]
            point1, point2 = ast.literal_eval(complement_info)
            b_chromosome = self.complement(bits_seq, point1, point2)
            b_chromosome = group_bits(b_chromosome)
            new_chromosome = bits_to_dna(b_chromosome, utils.two_bits_to_dna_base_table)

            new_population.append(new_chromosome)

        return new_population

    def dna_gdt(self, text, key):
        print("\nDNA-GDT is running...\n")

        rounds_no = int(get_pattern(no_rounds_del, key)[0])
        rounds = get_pattern(round_del, key)

        b_data = text
        print("Initial DNA sequence:", b_data)

        # run the algorithm "rounds_no" times
        while rounds_no > 0:
            round_info = rounds[rounds_no - 1]

            # create the chromosome population
            b_data = self.reshape(b_data, get_pattern(reshape_del, round_info))
            # print("Population data:", b_data)

            # apply mutation on population
            b_data = self.mutation(b_data, get_pattern(mutation_del, round_info))
            # print("Population data:", b_data)

            # apply crossover on population
            b_data = self.crossover(b_data, round_info)
            # print("Population data:", b_data)

            # decrypt data with key after reshaping it back to binary sequence and then convert it back to dna sequence
            # where decrypt = encrypt(encrypt(data, key), key) and encrypt => xor operation, because (a xor b) xor b = a
            encryption_key = get_pattern(key_del, key)[0]
            b_data = bits_to_dna(
                group_bits(
                    self.encrypt_key(dna_to_bits(self.reverse_reshape(b_data), utils.dna_base_to_two_bits_table),
                                     encryption_key)),
                utils.two_bits_to_dna_base_table)
            # print("Decrypted data:", b_data)

            rounds_no -= 1

        return bin2str(dna_to_bits(b_data, utils.dna_base_to_two_bits_table)).decode()

""" --------------------------------------------------------------------- """
""" --------------------------------------------------------------------- """


class DNA_Cryptography_Encryption(QMainWindow):
    def __init__(self):
        super(DNA_Cryptography_Encryption, self).__init__()

        loadUi("DNAEncryption.ui", self)

        self.update_files()

        self.dnaEncFunc_Obj = DNA_Encryption_Functions()
        self.dnaEncFunc_Obj.set_globals()

        self.dnaDcrFunc_Obj = DNA_Decryption_Functions()

        self.dna_key = str2bin(
            ''.join(random.SystemRandom().choice(string.ascii_letters + string.digits) for _ in range(16)))

        global decryption_key
        # append the key first
        decryption_key += key_del + self.dna_key + key_del

        utils.generate_pre_processing_tables()
        utils.generate_mutation_tables()

        self.key0 = DesKey(b"1234567812345678REAL_KEY")
        self.dna_encrypted_text = ''
        self.encrypted_msg_des = ''

        self.encryption_pushButton.clicked.connect(self.Encryption_Function)
        self.decrytion_pushButton.clicked.connect(self.Decryption_Function)
        self.exit_pushButton.clicked.connect(self.Exit_Function)
        self.algo_comboBox.activated.connect(self.update_files)

    @pyqtSlot()
    def update_files(self):

        # if os.path.isfile(key_filename):
        #     os.remove(key_filename)
        if os.path.isfile(encrypted_filename):
            os.remove(encrypted_filename)

    @pyqtSlot()
    def Encryption_Function(self):
        text = self.input_plainTextEdit.toPlainText()

        if len(text) != 0:

            print()
            print(text)
            print()

            start = time()
            self.dna_encrypted_text = self.dnaEncFunc_Obj.dna_get(text, self.dna_key)
            print("Final DNA sequence:", self.dna_encrypted_text)
            end = time()

            print("\nTotal execution time:", end - start)

            key_file = open(key_filename, "w")
            encrypted_file = open(encrypted_filename, "w")

            # save the encryption to a file to be used in the decryption process
            encrypted_file.write(self.dna_encrypted_text)

            self.dna_encrypted_plainTextEdit.setPlainText(str(self.dna_encrypted_text))

            # save key to a file to be read in the decryption process
            key_file.write(decryption_key)

            encrypted_file.close()
            key_file.close()

            algo_select = self.algo_comboBox.currentText()

            if algo_select == 'DES':

                message = self.dna_encrypted_text
                self.encrypted_msg_des = self.key0.encrypt(str.encode(message), padding=True)

                print("original string: ", message)
                print("encrypted string: ", self.encrypted_msg_des)

                self.algo_encrypted_plainTextEdit.setPlainText(str(self.encrypted_msg_des))


    @pyqtSlot()
    def Decryption_Function(self):
        algo_select = self.algo_comboBox.currentText()

        key_file = open(key_filename, "r")
        encrypted_file = open(encrypted_filename, "r")
        decrypted_file = open(decrypted_filename, "w")

        encrypted_text = encrypted_file.read()
        key = key_file.read()

        print("Encrypted text:", encrypted_text)
        print("Key:", key)

        # generate the encoding tables
        utils.generate_pre_processing_tables()
        utils.generate_mutation_tables()

        decrypted_algo_message = ''

        if algo_select == 'DES':

            decrypted_msg = self.key0.decrypt(self.encrypted_msg_des, padding=True)
            decrypted_msg = decrypted_msg.decode()
            print("decrypted string: ", decrypted_msg)

            decrypted_algo_message = decrypted_msg

            self.dna_decrypted_plainTextEdit.setPlainText(decrypted_msg)        

        print(decrypted_algo_message)

        # get the decrypted text
        start = time()
        decrypted_text = self.dnaDcrFunc_Obj.dna_gdt(decrypted_algo_message, key)
        print("Decrypted text:", decrypted_text)
        decrypted_file.write(decrypted_text)
        end = time()

        self.original_text_plainTextEdit.setPlainText(str(decrypted_text))

        encrypted_file.close()
        decrypted_file.close()
        key_file.close()

    def Exit_Function(self):
        sys.exit()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = DNA_Cryptography_Encryption()
    window.show()
    sys.exit(app.exec_())
